import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Kayser2005")

@time begin
    import SparseArrays
    import Base.Threads: @threads, threadid

    # -------------------------------------------------------------------
    # Run add https://github.com/josePereiro/Chemostat_Kayser2005.jl in 
    # the Julia Pkg REPL to install the package, then you must activate 
    # the package enviroment (see README)
    import Chemostat_Kayser2005
    const ChK = Chemostat_Kayser2005
    const iJR = ChK.iJR904
    const Kd = ChK.KayserData # experimental data
    const Bd = ChK.BegData    # cost data

    # -------------------------------------------------------------------
    # run add "https://github.com/josePereiro/Chemostat" in the 
    # julia Pkg REPL for installing the package
    import Chemostat
    import Chemostat.LP.MathProgBase

    const Ch = Chemostat
    const ChP = Ch.Plots
    const ChU = Ch.Utils
    const ChSS = Ch.SteadyState
    const ChLP = Ch.LP
    const ChEP = Ch.MaxEntEP
    const ChSU = Ch.SimulationUtils

    import UtilsJL
    const UJL = UtilsJL

    using Serialization

    # -------------------------------------------------------------------
    using Plots, FileIO
    import GR
    GR.inline("png")

end

## -------------------------------------------------------------------
INDEX = ChU.load_data(iJR.MAXENT_VARIANTS_INDEX_FILE; verbose = false);

## -------------------------------------------------------------------
const HOMO = :HOMO
const BOUNDED = :BOUNDED
const EXPECTED = :EXPECTED
const FBA = :FBA

## -------------------------------------------------------------------
fileid = "2.1"
function mysavefig(p, pname; params...) 
    fname = UJL.mysavefig(p, string(fileid, "_", pname), iJR.MODEL_FIGURES_DIR; params...)
    @info "Plotting" fname
end
myminmax(a::Vector) = (minimum(a), maximum(a))
CONC_IDERS = ["GLC", "AC", "NH4"]
FLX_IDERS = ["GLC", "CO2", "O2", "AC", "NH4"]

EXPS = 1:13 # experiments that have both concentration and flx data

exp_colors = let
    colors = Plots.distinguishable_colors(length(EXPS))
    Dict(exp => color for (exp, color) in zip(EXPS, colors))
end

ider_colors = let
    colors = Plots.distinguishable_colors(length(FLX_IDERS))
    Dict(met => color for (met, color) in zip(FLX_IDERS, colors))
end

method_colors = Dict(
    HOMO => :red,
    BOUNDED => :orange,
    EXPECTED => :blue,
)

## -------------------------------------------------------------------
# Collect
DAT = ChU.DictTree()
let 
    objider = iJR.KAYSER_BIOMASS_IDER
    exch_met_map = iJR.load_exch_met_map()
    Kd_mets_map = iJR.load_mets_map()

    for method in [HOMO, EXPECTED, BOUNDED]
        for exp in EXPS
            
            # !haskey(INDEX, method, :DFILE, exp) && continue
            datfile = INDEX[method, :DFILE, exp]
            dat = deserialize(datfile)
            
            global dat0 = dat
            
            model = dat[:model]
            objidx = ChU.rxnindex(model, objider)
            epouts = dat[:epouts]
            exp_beta = maximum(keys(epouts))
            epout = epouts[exp_beta]
            exp_xi = Kd.val(:xi, exp)

            println()
            @info("Doing", exp, method, length(dat[:epouts]), epout.iter); 
            # fbaout = dat[:fbaout]

            # Biomass
            # fba_biom = ChU.av(model, fbaout, objidx)
            ep_biom = ChU.av(model, epout, objidx)
            ep_std = sqrt(ChU.va(model, epout, objidx))
            Kd_biom = Kd.val("D", exp)
            
            # store
            # DAT[:fba  , :flx, objider, exp] = fba_biom
            DAT[method, :ep   , :flx, objider, exp] = ep_biom
            DAT[method, :eperr, :flx, objider, exp] = ep_std
            DAT[method, :Kd   , :flx, objider, exp] = Kd_biom
            DAT[:Kd   , :flx, objider, exp] = Kd_biom
            DAT[method, :fva  , :flx, objider, exp] = ChU.bounds(model, objider)
            
            # fuxes
            for Kd_met in FLX_IDERS

                    model_met = Kd_mets_map[Kd_met]
                    model_exch = exch_met_map[model_met]
                    model_exchi = ChU.rxnindex(model, model_exch)

                    ep_av = ChU.av(model, epout, model_exchi)
                    ep_std = sqrt(ChU.va(model, epout, model_exchi))
                    Kd_flx = Kd.val("u$Kd_met", exp)
                    
                    DAT[method, :Kd, :flx, Kd_met, exp] = Kd_flx
                    DAT[:Kd, :flx, Kd_met, exp] = Kd_flx
                    DAT[method, :ep, :flx, Kd_met, exp] = ep_av
                    DAT[method, :eperr, :flx, Kd_met, exp] = ep_std
                    
                    DAT[method, :fva , :flx, Kd_met, exp] = ChU.bounds(model, model_exch)

            end

            # mets
            for Kd_met in CONC_IDERS

                ep_std = DAT[method, :eperr, :flx, Kd_met, exp] 
                ep_av = DAT[method, :ep, :flx, Kd_met, exp]
                # conc (s = c + u*xi)
                c = Kd.val("c$Kd_met", exp, 0.0)
                # fba_conc = max(c + fba_av * exp_xi, 0.0)
                ep_conc = max(c + ep_av * exp_xi, 0.0)
                Kd_conc = Kd.val("s$Kd_met", exp)

                DAT[method, :Kd, :conc, Kd_met, exp] = Kd_conc
                DAT[:Kd, :conc, Kd_met, exp] = Kd_conc
                DAT[method, :ep, :conc, Kd_met, exp] = ep_conc
                DAT[method, :eperr, :conc, Kd_met, exp] = ep_std * exp_xi
            end

        end # for exp in EXPS
    
    end # for method

end

## -------------------------------------------------------------------
# beta vs stuff
let
    method = EXPECTED
    cGLC_plt = plot(;xlabel = "cGLC", ylabel = "beta")
    D_plt = plot(;xlabel = "D", ylabel = "beta")
    for exp in EXPS 
        datfile = INDEX[method, :DFILE, exp]
        dat = deserialize(datfile)
        beta = maximum(keys(dat[:epouts]))

        params = (;label = "", color = exp_colors[exp], 
            alpha = 0.7, ms = 7
        )
        cGLC = Kd.val("cGLC", exp)
        D = Kd.val("D", exp)
        scatter!(cGLC_plt, [cGLC], [beta]; params...)
        scatter!(D_plt, [D], [beta]; params...)
    end
    mysavefig([cGLC_plt, D_plt], "beta_vs_stuff"; method)
end

## -------------------------------------------------------------------
# EP biomass corr
let
    objider = iJR.KAYSER_BIOMASS_IDER
    ps = Plots.Plot[]
    for method in [HOMO, EXPECTED, BOUNDED]
        p = plot(title = string(iJR.PROJ_IDER, " method: ", method), 
            xlabel = "model biom", ylabel = "exp biom")
        ep_vals = DAT[method, :ep, :flx, objider, EXPS]
        eperr_vals = DAT[method, :eperr, :flx, objider, EXPS]
        Kd_vals = DAT[method, :Kd, :flx, objider, EXPS]
        color = [exp_colors[exp] for exp in EXPS]
        m, M = myminmax([Kd_vals; ep_vals])
        margin = abs(M - m) * 0.1
        scatter!(p, ep_vals, Kd_vals; 
            xerr = eperr_vals,
            label = "", color,
            alpha = 0.7, ms = 7,
            xlim = [m - margin, M + margin],
            ylim = [m - margin, M + margin],
        )
        push!(ps, p)
    end
    layout = (1, length(ps))
    mysavefig(ps, "obj_val_ep_corr"; layout)
end

## -------------------------------------------------------------------
# EXPECTED flux vs beta
let
    objider = iJR.KAYSER_BIOMASS_IDER
    method = EXPECTED
    p = plot(title = iJR.PROJ_IDER, xlabel = "beta", ylabel = "biom")
    for exp in EXPS 
        datfile = INDEX[method, :DFILE, exp]
        dat = deserialize(datfile)
        model = dat[:model]
        objidx = ChU.rxnindex(model, objider)
        epouts = dat[:epouts]
        exp_beta = maximum(keys(epouts))
        exp_xi = Kd.val("xi", exp)
        scatter!(p, [exp_beta], [Kd.val("D", exp)], ms = 12, color = :white, label = "")

        betas = collect(keys(epouts)) |> sort
        bioms = [ChU.av(model, epouts[beta], objidx) for beta in betas]
        scatter!(p, betas, bioms, label = "", color = :black, alpha = 0.2)

    end
    mysavefig(p, "obj_val_vs_beta"; method)
end

## -------------------------------------------------------------------
# total correlations
let
    for (dat_prefix, iders, zoom_lim) in [(:flx, FLX_IDERS, [-2.5, 2.5]), 
                                            (:conc, CONC_IDERS, [0.0, 100.0])]

        ps = Plots.Plot[]
        for method in [HOMO, EXPECTED, BOUNDED]                                            
            ep_vals = DAT[method, :ep, dat_prefix, iders, EXPS]
            ep_errs = DAT[method, :eperr, dat_prefix, iders, EXPS]
            Kd_vals = DAT[method, :Kd, dat_prefix, iders, EXPS]
            color = [ider_colors[ider] for ider in iders, exp in EXPS]
            m, M = myminmax([ep_vals; Kd_vals])

            scatter_params = (;label = "", color, ms = 7, alpha = 0.7)
            # ep corr
            p1 = plot(title = "$(iJR.PROJ_IDER) (EP) $method", 
                ylabel = "model $(dat_prefix)",
                xlabel = "exp $(dat_prefix)", 
            )
            scatter!(p1, Kd_vals, ep_vals; xerr = ep_errs, scatter_params...)
            plot!(p1, [m,M], [m,M]; ls = :dash, color = :black, label = "")
            push!(ps, deepcopy(p1))

        end

        # fba corr
        # # fba_vals = DAT[:fba, dat_prefix, iders, EXPS]
        # p2 = plot(title = "$(iJR.PROJ_IDER) (FBA)", 
        #     ylabel = "exp $(dat_prefix)",
        #     xlabel = "model $(dat_prefix)", 
        # )
        # scatter!(p2, fba_vals, Kd_vals; scatter_params...)
        # plot!(p2, [m,M], [m,M]; ls = :dash, color = :black, label = "")

        # mysavefig([p1, p2], pname)

        # # zoom
        # pname = string(dat_prefix, "_tot_corr_zoom")
        # ps = plot.([p1, p2]; xlim = zoom_lim, ylim = zoom_lim)
        # mysavefig(ps, pname)

        layout = (1, length(ps))
        pname = string(dat_prefix, "_tot_corr")
        mysavefig(ps, pname; layout)
    end

end

## -------------------------------------------------------------------
# fva bounds
let
   
    ps = Plots.Plot[]
    for ider = FLX_IDERS
        p = plot(title = ider, xlabel = "replica", ylabel = "flx")
        xticks =  (EXPS, string.(EXPS))
        
        Kd_vals = DAT[:Kd, :flx, ider, EXPS]
        plot!(p, EXPS, Kd_vals; 
            label = "exp", color = :black, alpha = 0.8, lw = 3, xticks)

        for method in [HOMO, EXPECTED, BOUNDED]             
            color = method_colors[method]    
            
            ep_vals = DAT[method, :ep, :flx, ider, EXPS]
            plot!(p, EXPS, ep_vals; 
                label = string(method), color, alpha = 0.5, lw = 5, ls = :dash, xticks)
            
            # fba_vals = DAT[method, :fba, :flx, ider, EXPS]
            # plot!(p, EXPS, fba_vals; 
            #     label = "fba", color = :red, alpha = 0.5, lw = 5, ls = :dash, xticks)
            
            fva_ranges = DAT[method, :fva, :flx, ider, EXPS]
            plot!(p, EXPS, last.(fva_ranges);  
                label = "", color, alpha = 0.8, ls = :dot, lw = 3, xticks)
            plot!(p, EXPS, first.(fva_ranges); 
                label = "", color, alpha = 0.8, ls = :dot, lw = 3, xticks)
        end
        push!(ps, p)
    end
    pname = string("bound_study")
    mysavefig(ps, pname)
    
end


## -------------------------------------------------------------------
# let
#     Kd_met, exp = "GLC", 1
#     datfile = INDEX[method, :DFILE, exp]
#     dat = deserialize(datfile)
#     # dat[[:ep, :fba, :Kd], :flx, Kd_met, exp]
#     dat[:ep, :flx, Kd_met, exp]
# end

## -------------------------------------------------------------------
# marginal distributions
let 
    objider = iJR.BIOMASS_IDER
    size = [300, 250]
    # method = EXPECTED

    # Iders
    model_iders, Kd_iders = [iJR.BIOMASS_IDER], ["D"]
    for Kd_met in CONC_IDERS
        model_met = iJR.Kd_mets_map[Kd_met]
        model_exch = iJR.exch_met_map[model_met]
        push!(model_iders, model_exch)
        push!(Kd_iders, string("u", Kd_met))
    end
    
    for (model_ider, Kd_ider) in zip(model_iders, Kd_iders)
        ps = Plots.Plot[]
        ps2 = Plots.Plot[]
        for exp in EXPS
            p = plot(title = string(Kd_ider, " exp: ", exp))
            p2 = plot(title = string(Kd_ider, " exp: ", exp))
            margin, m, M = -Inf, Inf, -Inf
            Kd_av = Kd.val(Kd_ider, exp)
            
            # EP
            for method in [BOUNDED, EXPECTED, HOMO]
                color = method_colors[method]    

                datfile = INDEX[method, :DFILE, exp]
                dat = deserialize(datfile)
                model = dat[:model]
                objidx = ChU.rxnindex(model, objider)
                exp_beta = dat[:exp_beta]
                epouts = dat[:epouts]
                epout = epouts[exp_beta]
                ep_av = ChU.av(model, epout, model_ider)
                ep_va = sqrt(ChU.va(model, epout, model_ider))
                # fbaout = dat[:fbaout]
                        
                # ChP.plot_marginal!(p, model, [epout, fbaout], model_exch; legend = false)
                ChP.plot_marginal!(p, model, [epout], model_ider; 
                    legend = false, color, alpha = 0.6, lw = 5)
                
                m = minimum([m, ep_av, Kd_av])
                M = maximum([M, ep_av, Kd_av])
                margin = maximum([margin, 3 * ep_va])

                if method == EXPECTED
                    for (beta, epout) in sort(epouts; by = first)
                        ep_av = ChU.av(model, epout, model_ider)
                        ep_va = sqrt(ChU.va(model, epout, model_ider))
                        
                        alpha = 0.15
                        color = method_colors[method]
                        ChP.plot_marginal!(p2, model, epout, model_ider; 
                            legend = false, color, alpha, lw = 1)

                        if beta == exp_beta
                            # @info "At" Kd_ider exp
                            ChP.plot_marginal!(p2, model, epout, model_ider; 
                                legend = false, color, 
                                alpha = 1.0, lw = 3
                            )
                            break
                        end
                    end
                    push!(ps2, p2)
                end

            end
            # Experimental
            vline!(p, [Kd_av]; label = "", lw = 6, color = :black, alpha = 0.3)
            vline!(p2, [Kd_av]; label = "", lw = 6, color = :black, alpha = 0.3)
            
            plot!(p; xlim = [m - margin, M + margin], size)
            plot!(p2; xlim = [m - margin, M + margin], size)
            push!(ps, p)
        end

        for k in [:xi, :D, :cGLC]
            p = plot(;title = Kd_ider, size)
            xticks =  (EXPS, string.(EXPS))
            p = bar!(p, EXPS, Kd.val(k); title = k, label = "", xticks)
            push!(ps, p)
            push!(ps2, p)
        end

        pname = string(Kd_ider, "_marginals")
        mysavefig(ps, pname)
        pname = string(Kd_ider, "_marginals_vs_beta")
        mysavefig(ps2, pname)
    end

end 

## -------------------------------------------------------------------
# leyends
# TODO fix this...
let
    for (title, colors) in [
            ("exp", exp_colors), 
            ("iders", ider_colors),
            ("method", method_colors)
        ]
    p = plot(; framestyle = :none)
        scolors = sort(collect(colors); by = (p) -> string(first(p)))
        for (id, color) in scolors
            scatter!(p, [0], [0];
                thickness_scaling = 1,
                color, ms = 8, label = string(id),
                legendfontsize=10, 
                # size = [300, 900],
                # legend = :left
            )
        end
        mysavefig(p, "$(title)_color_legend")
    end
end