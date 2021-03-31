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
    const ChU = Ch.Utils
    const ChSS = Ch.SteadyState
    const ChLP = Ch.LP
    const ChEP = Ch.MaxEntEP
    const ChSU = Ch.SimulationUtils

    import ChemostatPlots
    const ChP = ChemostatPlots
    
    import UtilsJL
    const UJL = UtilsJL

    using ProgressMeter
    using Serialization
    using Statistics

    # -------------------------------------------------------------------
    using Plots, FileIO
    import GR
    GR.inline("png")
    using Base.Threads

end

## -------------------------------------------------------------------
INDEX = ChU.load_data(iJR.MAXENT_VARIANTS_INDEX_FILE; verbose = false);

# -------------------------------------------------------------------
const ME_Z_OPEN_G_OPEN         = :ME_Z_OPEN_G_OPEN
const ME_MAX_POL               = :ME_MAX_POL
const ME_Z_EXPECTED_G_MOVING   = :ME_Z_EXPECTED_G_MOVING
const ME_Z_EXPECTED_G_BOUNDED  = :ME_Z_EXPECTED_G_BOUNDED
const ME_Z_FIXXED_G_BOUNDED    = :ME_Z_FIXXED_G_BOUNDED

ALL_METHODS = [
    # ME_Z_OPEN_G_OPEN, 
    ME_MAX_POL,
    # ME_Z_FIXXED_G_BOUNDED, 
    # ME_Z_EXPECTED_G_BOUNDED, 
    # ME_Z_EXPECTED_G_MOVING
]

# -------------------------------------------------------------------
fileid = "2.1"
function mysavefig(p, pname; params...)
    pname = string(fileid, "_", pname)
    fname = UJL.mysavefig(p, pname, iJR.MODEL_FIGURES_DIR; params...)
    @info("Plotting", fname)
end
myminmax(a::Vector) = (minimum(a), maximum(a))
# CONC_IDERS = ["GLC", "AC", "NH4"]
FLX_IDERS = ["GLC", "CO2", "O2", "AC", "NH4"]
# FLX_IDERS = ["GLC", "AC", "NH4"]

Kd_rxns_map = iJR.load_Kd_rxns_map()
Kd_mets_map = iJR.load_Kd_mets_map()

## -------------------------------------------------------------------
# Collect
DAT = ChU.DictTree()
let 

    WLOCK = ReentrantLock()

    # CACHE
    DATfile = joinpath(iJR.MODEL_PROCESSED_DATA_DIR, "2.1_DAT.jls")
    # if isfile(DATfile) 
    #     global DAT = deserialize(DATfile) 
    #     @info("DAT CACHE LOADED")
    #     return
    # end
    DAT[:EXPS] = []

    objider = iJR.KAYSER_BIOMASS_IDER
    DAT[:CONC_IDERS] = CONC_IDERS
    DAT[:FLX_IDERS] = FLX_IDERS

    # Find exps
    for exp in Kd.EXPS
        ok = false
        for method in ALL_METHODS
            ok = haskey(INDEX, method, :DFILE, exp) &&
                INDEX[method, :DFILE, exp] != :unfeasible
            !ok && break
        end
        !ok && continue
        push!(DAT[:EXPS], exp)
    end
    max_model = iJR.load_model("max_model"; uncompress = false)

    # Feed jobs
    nths = nthreads()
    Ch = Channel(nths) do ch
        for exp in DAT[:EXPS], method in ALL_METHODS
            put!(ch, (exp, method))
        end
    end


    @threads for thid in 1:nths
        thid = threadid()
        thid == 1 && nths > 1 && continue
        for (exp, method) in Ch
                
            !haskey(INDEX, method, :DFILE, exp) && continue
            datfile = INDEX[method, :DFILE, exp]
            datfile == :unfeasible && continue
            dat = deserialize(datfile)
            
            model = dat[:model]
            objidx = ChU.rxnindex(model, objider)
            epouts = dat[:epouts]
            exp_beta = maximum(keys(epouts)) # dat[:exp_beta]
            epout = epouts[exp_beta]
            exp_xi = Kd.val(:xi, exp)
            fva_model = iJR.load_model("fva_models", exp; uncompress = false)
            
            lock(WLOCK) do
                @info("Doing", 
                    exp, method, 
                    length(dat[:epouts]), 
                    epout.iter, thid
                ); println()
            end

            # Biomass
            ep_biom = ChU.av(model, epout, objidx)
            ep_std = sqrt(ChU.va(model, epout, objidx))
            Kd_biom = Kd.val("D", exp)
            max_lb, max_ub = ChU.bounds(max_model, objidx)
            fva_lb, fva_ub = ChU.bounds(fva_model, objidx)
            lb = max(max_lb, fva_lb)
            ub = min(max_ub, fva_ub)
            
            # store
            lock(WLOCK) do
                DAT[method, :ep   , :flx, objider, exp] = ep_biom
                DAT[method, :eperr, :flx, objider, exp] = ep_std
                DAT[method, :Kd   , :flx, objider, exp] = Kd_biom
                DAT[:Kd   , :flx, objider, exp] = Kd_biom
                DAT[method, :bounds, :flx, objider, exp] = (lb, ub)
            end

            # fluxes
            for Kd_met in FLX_IDERS

                    model_met = Kd_mets_map[Kd_met]
                    model_exch = Kd_rxns_map[Kd_met]
                    model_exchi = ChU.rxnindex(model, model_exch)

                    ep_av = ChU.av(model, epout, model_exchi)
                    ep_std = sqrt(ChU.va(model, epout, model_exchi))
                    Kd_flx = Kd.val("u$Kd_met", exp)
                    proj = ChLP.projection2D(model, objider, model_exchi; l = 50)
                    
                    max_lb, max_ub = ChU.bounds(max_model, Kd_rxns_map[Kd_met])
                    fva_lb, fva_ub = ChU.bounds(fva_model, Kd_rxns_map[Kd_met])
                    lb = max(max_lb, fva_lb)
                    ub = min(max_ub, fva_ub)
                            
                    lock(WLOCK) do
                        DAT[method, :Kd, :flx, Kd_met, exp] = Kd_flx
                        DAT[:Kd, :flx, Kd_met, exp] = Kd_flx
                        DAT[method, :ep, :proj, Kd_met, exp] = proj
                        DAT[method, :ep, :flx, Kd_met, exp] = ep_av
                        DAT[method, :eperr, :flx, Kd_met, exp] = ep_std
                        DAT[method, :bounds, :flx, Kd_met, exp] = (lb, ub)
                    end

            end

            # mets
            for Kd_met in CONC_IDERS

                ep_std = DAT[method, :eperr, :flx, Kd_met, exp] 
                ep_av = DAT[method, :ep, :flx, Kd_met, exp]
                # conc (s = c + u*xi)
                c = Kd.val("c$Kd_met", exp, 0.0)
                ep_conc = max(c + ep_av * exp_xi, 0.0)
                Kd_conc = Kd.val("s$Kd_met", exp)

                lock(WLOCK) do
                    DAT[method, :Kd, :conc, Kd_met, exp] = Kd_conc
                    DAT[:Kd, :conc, Kd_met, exp] = Kd_conc
                    DAT[method, :ep, :conc, Kd_met, exp] = ep_conc
                    DAT[method, :eperr, :conc, Kd_met, exp] = ep_std * exp_xi
                end
            end

        end # for (exp, method)
    end # for thid

    # saving
    DAT[:EXPS] |> unique! |> sort!
    serialize(DATfile, DAT)
end
EXPS = DAT[:EXPS]

exp_colors = let
    colors = Plots.distinguishable_colors(length(EXPS))
    Dict(exp => color for (exp, color) in zip(EXPS, colors))
end

ider_colors = Dict(
    "GLC" => :red, "CO2" => :yellow,
    "O2" => :blue, "AC" => :orange, 
    "NH4" => :green, "D" => :black,
)

method_colors = Dict(
    ME_Z_OPEN_G_OPEN => :red,
    ME_MAX_POL => :blue,
    ME_Z_EXPECTED_G_BOUNDED => :orange,
    ME_Z_EXPECTED_G_MOVING => :purple,
    ME_Z_FIXXED_G_BOUNDED => :green,
)

# -------------------------------------------------------------------
# Inter project comunication
let
    CORR_DAT = isfile(iJR.CORR_DAT_FILE) ? ChU.load_data(iJR.CORR_DAT_FILE) : Dict()
    CORR_DAT[:MAXENT_EP] = DAT
    ChU.save_data(iJR.CORR_DAT_FILE, CORR_DAT)
end

##  ----------------------------------------------------------------------------
# pol vox volumen
let
    model0 = iJR.load_model("max_model")

    bins = 50
    Kd_rxns_map = iJR.load_Kd_rxns_map()
    offsetf = 1.1
    
    exglcidx = ChU.rxnindex(model0, iJR.GLC_EX_IDER)
    biomidx = ChU.rxnindex(model0, iJR.KAYSER_BIOMASS_IDER)
    exglcL, exglcU = ChU.bounds(model0, exglcidx)
    
    maxD = maximum(Kd.val(:D)) * offsetf 
    max_cgD_X = -maximum(Kd.ciD_X(:GLC)) * offsetf
    
    Ds = range(0.01, maxD; length = bins)
    cgD_Xs = range(max_cgD_X, exglcU; length = bins)

    box_vols = zeros(bins, bins) 
    
    Kd_ids = ["GLC", "AC", "NH4"]
    model_ids = [Kd_rxns_map[Kd_id] for Kd_id in Kd_ids]
    model_idxs = [ChU.rxnindex(model0, model_id) for model_id in model_ids]

    # feeding task
    nths = nthreads()
    Ch = Channel(nths) do Ch_
        @showprogress for (Di, D) in enumerate(Ds)
            for (cgD_Xi, cgD_X) in enumerate(cgD_Xs)
                put!(Ch_, (Di, D, cgD_Xi, cgD_X))
            end
        end 
    end

    # compute volume map
    @threads for _ in 1:nths
        thid = threadid()
        thid == 1 && nths > 1 && continue
        model = deepcopy(model0)
        for (Di, D, cgD_Xi, cgD_X) in Ch
            
            # Reduce Pol
            ChU.lb!(model, exglcidx, cgD_X)
            ChU.bounds!(model, biomidx, D, D)
            try
                L, U = ChLP.fva(model, model_idxs)
                vol = prod(abs.(L .- U))
                box_vols[Di, cgD_Xi] = max(0.0, log10(vol + 1e-50))
            catch
                box_vols[Di, cgD_Xi] = NaN
            end
        end
    end

    # vol map
    p = heatmap(Ds, cgD_Xs, box_vols'; 
        title = "Polytope Box volume", label = "", 
        xlabel = "D", ylabel = "cgD/ X"
    )

    # exp vals
    exp_Ds = [Kd.val(:D, exp) for exp in EXPS]
    exp_cgD_Xs = [-Kd.ciD_X(:GLC, exp) for exp in EXPS]
    scatter!(p, exp_Ds, exp_cgD_Xs;
        label = "exp", color = :white, 
        m = 8
    )

    mysavefig(p, "pol_box_volume"; bins)
end

## -------------------------------------------------------------------
# Var
let
    method = ME_MAX_POL

    var_avs = Dict()
    var_stds = Dict()
    for exp in EXPS
        datfile = INDEX[method, :DFILE, exp]
        dat = deserialize(datfile)
        
        epouts = dat[:epouts]
        exp_beta = maximum(keys(epouts))
        epout = epouts[exp_beta]
        
        vas = ChU.va(epout)
        var_avs[exp] = mean(vas)
        var_stds[exp] = 0.0 #std(vas)
    end

    # sGLC
    ids = [:D, :sAC, :sGLC, :sNH4, :Xv, :uAC, :uGLC, :uNH4]
    for id in ids
        p = scatter(;xlabel = string(id), ylabel = "avs log var")
        xs = [Kd.val(id, exp) for exp in EXPS]
        ys = [var_avs[exp] for exp in EXPS]
        errs = [var_stds[exp] for exp in EXPS]
        sids = sortperm(xs)
        ys, errs = ys[sids], errs[sids]
        scatter!(p, xs, ys; yerr = errs, label = "", color = :black, m = 8)
        plot!(p, xs, ys; label = "", color = :black, lw = 3, alpha = 0.5)
        mysavefig(p, string("av_log_var_vs_", id))
    end
end

## -------------------------------------------------------------------
# MSE per method
let

    p = plot(;xlabel = "experiment", ylabel = "MSE")
    for method in ALL_METHODS
        MSEs = []
        for exp in EXPS

            sum = 0.0
            N = 0
            glc_flx = DAT[method, :Kd, :flx, "GLC", exp]
            for ider in FLX_IDERS
                model_val = DAT[method, :ep, :flx, ider, exp]
                exp_val = DAT[method, :Kd, :flx, ider, exp]
                sum += (model_val/glc_flx - exp_val/glc_flx)^2
                N += 1
            end
            push!(MSEs, sum / N)
        end
        scatter!(p, EXPS, MSEs; color = method_colors[method],
            label = string(method), m = 8, alpha = 0.8, 
            legend = :topleft
        )
        plot!(p, EXPS, MSEs; color = method_colors[method],
            label = "", ls = :dash, alpha = 0.8
        )
    end
    mysavefig(p, "MSE_per_method")
end

## -------------------------------------------------------------------
# MSE per ider
let
    p = plot(;xlabel = "experiment", ylabel = "MSE")
    for method in ALL_METHODS
        MSEs = []

        for ider in FLX_IDERS
            sum = 0.0
            N = 0
            for exp in EXPS
                glc_flx = DAT[method, :Kd, :flx, "GLC", exp]
                model_val = DAT[method, :ep, :flx, ider, exp]
                exp_val = DAT[method, :Kd, :flx, ider, exp]
                sum += (model_val/glc_flx - exp_val/glc_flx)^2
                N += 1
            end
            push!(MSEs, sum / N)
        end

        scatter!(p, FLX_IDERS, MSEs; color = method_colors[method],
            label = string(method), m = 8, alpha = 0.8, 
            legend = :topleft
        )
        plot!(p, FLX_IDERS, MSEs; color = method_colors[method],
            label = "", ls = :dash, alpha = 0.8
        )
    end
    mysavefig(p, "MSE_per_ider")
end

# ## -------------------------------------------------------------------
# # MSE per beta
# let
#     method = ME_Z_EXPECTED_G_BOUNDED

#     ps = Plots.Plot[]
#     for exp in EXPS
#         p = plot(;title = string("exp: ", exp), 
#             xlabel = "beta", ylabel = "MSE"
#         )

#         datfile = INDEX[method, :DFILE, exp]
#         dat = deserialize(datfile)
#         epouts = dat[:epouts]
#         betas = epouts |> keys |> collect |> sort
#         exp_beta = maximum(keys(epouts)) # dat[:exp_beta]
#         model = dat[:model]
        
#         MSEs = []
#         for beta in betas
#             epout = epouts[beta]

#             sum = 0.0
#             N = 0

#             glc_flx = Kd.uval(:GLC, exp)
#             for ider in FLX_IDERS

#                 model_met = Kd_mets_map[ider]
#                 model_exch = Kd_rxns_map[model_met]
#                 model_exchi = ChU.rxnindex(model, model_exch)

#                 model_flx = ChU.av(model, epout, model_exchi)
#                 exp_flx = Kd.uval(ider, exp)

#                 sum += (model_flx/glc_flx - exp_flx/glc_flx)^2
#                 N += 1
#             end
            
#             push!(MSEs, sum / N)
#         end

#         scatter!(p, betas, MSEs; color = :black,
#             label = "", m = 8, alpha = 0.8
#         )
#         plot!(p, betas, MSEs; color = :black,
#             label = "", ls = :dash, alpha = 0.8
#         )
#         vline!(p, [exp_beta]; color = :black, 
#             label = "", ls = :dot, lw = 3, alpha = 0.9
#         )
#         push!(ps, p)
#     end
#     mysavefig(ps, "MSE_vs_beta")
# end

## -------------------------------------------------------------------
# proj 2D
let
    method = ME_MAX_POL
    biom_ider = iJR.KAYSER_BIOMASS_IDER
    exO2_ider = Kd_rxns_map["O2"]

    ps_pool = Dict()
    for exp in EXPS

        for Kd_ider in FLX_IDERS

            # 2D Projection
            p = plot(;title = string("Kayser2005 exp:", exp), 
                xlabel = string(biom_ider), ylabel = string(Kd_ider),
                legend = :left
            )
            proj = DAT[method, :ep, :proj, Kd_ider, exp]
            ChP.plot_projection2D!(p, proj; l = 50)

            # bounds
            lb, ub =DAT[method, :bounds, :flx, Kd_ider, exp]
            hline!(p, [lb]; lw = 3, 
                label = "fva lb",
                color = :blue, ls = :solid
            )
            hline!(p, [ub]; lw = 3,
                label = "fva ub", 
                color = :red, ls = :solid
            )

            # EXPERIMENTAL FLXS
            exp_biom = DAT[method, :Kd, :flx, biom_ider, exp]
            exp_exch = DAT[method, :Kd, :flx, Kd_ider, exp]
            scatter!(p, [exp_biom], [exp_exch]; 
                m = 8, color = :red, label = "exp"
            )
            
            # MAXENT FLXS
            ep_biom = DAT[method, :ep, :flx, biom_ider, exp]
            ep_biom_err = DAT[method, :eperr, :flx, biom_ider, exp]
            ep_exch = DAT[method, :ep, :flx, Kd_ider, exp]
            ep_exch_err = DAT[method, :eperr, :flx, Kd_ider, exp]
            scatter!(p, [ep_biom], [ep_exch]; 
                xerr = [ep_biom_err], yerr = [ep_exch_err],
                m = 8, color = :blue, label = "maxent"
            )

            # mysavefig(p, "polytope"; Kd_ider, exp, method)
            ps_pool[(exp, Kd_ider)] = deepcopy(p)
        end
    end

    # collect 
    for exp in EXPS
        ps = Plots.Plot[ps_pool[(exp, Kd_ider)] for Kd_ider in FLX_IDERS]
        mysavefig(ps, "polytope"; exp, method)
    end

    for Kd_ider in FLX_IDERS
        ps = Plots.Plot[ps_pool[(exp, Kd_ider)] for exp in EXPS]
        mysavefig(ps, "polytope"; Kd_ider, method)
    end
end

## -------------------------------------------------------------------
# # beta vs stuff
# let
#     method = ME_Z_EXPECTED_G_MOVING
#     cGLC_plt = plot(;xlabel = "cGLC", ylabel = "beta")
#     D_plt = plot(;xlabel = "D", ylabel = "beta")
#     for exp in EXPS 
#         datfile = INDEX[method, :DFILE, exp]
#         dat = deserialize(datfile)
#         beta = maximum(keys(dat[:epouts]))

#         params = (;label = "", color = exp_colors[exp], 
#             alpha = 0.7, ms = 7
#         )
#         cGLC = Kd.val("cGLC", exp)
#         D = Kd.val("D", exp)
#         scatter!(cGLC_plt, [cGLC], [beta]; params...)
#         scatter!(D_plt, [D], [beta]; params...)
#     end
#     mysavefig([cGLC_plt, D_plt], "beta_vs_stuff"; method)
# end

## -------------------------------------------------------------------
# EP biomass corr
let
    objider = iJR.KAYSER_BIOMASS_IDER
    ps = Plots.Plot[]
    for method in ALL_METHODS
        p = plot(title = string(iJR.PROJ_IDER, " method: ", method), 
            xlabel = "exp biom", ylabel = "model biom")
        ep_vals = DAT[method, :ep, :flx, objider, EXPS]
        eperr_vals = DAT[method, :eperr, :flx, objider, EXPS]
        Kd_vals = DAT[method, :Kd, :flx, objider, EXPS]
        color = [exp_colors[exp] for exp in EXPS]
        m, M = myminmax([Kd_vals; ep_vals])
        margin = abs(M - m) * 0.1
        scatter!(p, Kd_vals, ep_vals; 
            yerr = eperr_vals,
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
# # flux vs beta
# let
#     objider = iJR.KAYSER_BIOMASS_IDER
#     method = ME_Z_EXPECTED_G_MOVING
#     p = plot(title = iJR.PROJ_IDER, xlabel = "beta", ylabel = "biom")
#     for exp in EXPS 
#         datfile = INDEX[method, :DFILE, exp]
#         dat = deserialize(datfile)
#         model = dat[:model]
#         objidx = ChU.rxnindex(model, objider)
#         epouts = dat[:epouts]
#         exp_beta = maximum(keys(epouts))
#         exp_xi = Kd.val("xi", exp)
#         scatter!(p, [exp_beta], [Kd.val("D", exp)], ms = 12, color = :white, label = "")

#         betas = collect(keys(epouts)) |> sort
#         bioms = [ChU.av(model, epouts[beta], objidx) for beta in betas]
#         scatter!(p, betas, bioms, label = "", color = :black, alpha = 0.2)

#     end
#     mysavefig(p, "obj_val_vs_beta"; method)
# end

## -------------------------------------------------------------------
# correlations
let

    tot_ps = Plots.Plot[]
    for method in ALL_METHODS        
        # total corr
        let            
            ep_vals = DAT[method, :ep, :flx, FLX_IDERS, EXPS] .|> abs
            ep_errs = DAT[method, :eperr, :flx, FLX_IDERS, EXPS] .|> abs
            Kd_vals = DAT[method, :Kd, :flx, FLX_IDERS, EXPS] .|> abs
            
            color = [ider_colors[ider] for ider in FLX_IDERS, exp in EXPS]
            scatter_params = (;label = "", color, ms = 7, alpha = 0.7)
            # ep corr
            p = plot(title = "$(iJR.PROJ_IDER) (EP) $method", 
                ylabel = "model abs flx",
                xlabel = "exp abs flx", 
            )
            scatter!(p, ep_vals, Kd_vals; yerr = ep_errs, scatter_params...)
            all_vals = [ep_vals; Kd_vals] |> sort!
            plot!(p, all_vals, all_vals; ls = :dash, color = :black, label = "")
            push!(tot_ps, deepcopy(p))
        end

        # per ider
        let       
            for ider in FLX_IDERS
                ep_vals = DAT[method, :ep, :flx, ider, EXPS] .|> abs
                ep_errs = DAT[method, :eperr, :flx, ider, EXPS] .|> abs
                Kd_vals = DAT[method, :Kd, :flx, ider, EXPS] .|> abs
                
                color = ider_colors[ider]
                scatter_params = (;label = "", color, ms = 7, alpha = 0.7)
                # ep corr
                p = plot(title = "$(iJR.PROJ_IDER) (EP) $method", 
                    ylabel = "model abs flx",
                    xlabel = "exp abs flx", 
                )
                scatter!(p, ep_vals, Kd_vals; yerr = ep_errs, scatter_params...)
                bounds = DAT[method, :bounds, :flx, ider, EXPS]
                lb, ub = minimum(first.(bounds)), maximum(last.(bounds))
                plot!(p, abs.([lb, ub]), abs.([lb, ub]); 
                    ls = :dash, color = :black, label = "", 
                    # xlim = [lb, ub], ylim = [lb, ub]
                )
                mysavefig(p, "corr"; ider, method)
            end
        end
    end

    layout = (1, length(tot_ps))
    mysavefig(tot_ps, "flx_tot_corr"; layout)

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

        for method in ALL_METHODS
            color = method_colors[method]    

            fva_ranges = DAT[method, :bounds, :flx, ider, EXPS]
            for (exp, (lb, ub)) in zip(EXPS, fva_ranges)
                plot!(p, [exp, exp], [lb, ub];  
                    label = "", color, alpha = 0.1, ls = :solid, 
                    lw = 50
                )
            end
            
            ep_vals = DAT[method, :ep, :flx, ider, EXPS]
            plot!(p, EXPS, ep_vals; 
                label = string(method), color, alpha = 0.5, 
                lw = 5, ls = :dash, xticks
            )
            
        end
        push!(ps, p)
    end
    pname = string("bound_study")
    mysavefig(ps, pname)
    
end

## -------------------------------------------------------------------
# marginal distributions
let 

    method2 = ME_MAX_POL

    objider = iJR.KAYSER_BIOMASS_IDER
    size = [300, 250]

    # Iders
    model_iders, Kd_iders = [objider], ["D"]
    for Kd_met in CONC_IDERS
        model_met = Kd_mets_map[Kd_met]
        model_exch = Kd_rxns_map[Kd_met]
        push!(model_iders, model_exch)
        push!(Kd_iders, string("u", Kd_met))
    end
    
    for (model_ider, Kd_ider) in zip(model_iders, Kd_iders)
        ps = Plots.Plot[]
        ps_bs = Plots.Plot[]
        for exp in EXPS
            p = plot(title = string(Kd_ider, " exp: ", exp))
            p_bs = plot(title = string(Kd_ider, " exp: ", exp))
            margin, m, M = -Inf, Inf, -Inf
            Kd_av = Kd.val(Kd_ider, exp)
            
            # EP
            for method in ALL_METHODS
                color = method_colors[method]    

                datfile = INDEX[method, :DFILE, exp]
                dat = deserialize(datfile)
                model = dat[:model]
                objidx = ChU.rxnindex(model, objider)
                epouts = dat[:epouts]
                exp_beta = maximum(keys(epouts))
                epout = epouts[exp_beta]
                ep_av = ChU.av(model, epout, model_ider)
                ep_va = sqrt(ChU.va(model, epout, model_ider))
                        
                ChP.plot_marginal!(p, model, [epout], model_ider; 
                    legend = false, color, alpha = 0.6, lw = 5)
                
                m = minimum([m, ep_av, Kd_av])
                M = maximum([M, ep_av, Kd_av])
                margin = maximum([margin, 3 * ep_va])

                if method == method2
                    for (beta, epout) in sort(epouts; by = first)
                        ep_av = ChU.av(model, epout, model_ider)
                        ep_va = sqrt(ChU.va(model, epout, model_ider))
                        
                        alpha = 0.15
                        color = method_colors[method]
                        ChP.plot_marginal!(p_bs, model, epout, model_ider; 
                            legend = false, color, alpha, lw = 1)

                        if beta == exp_beta
                            ChP.plot_marginal!(p_bs, model, epout, model_ider; 
                                legend = false, color, 
                                alpha = 1.0, lw = 3
                            )
                            break
                        end
                    end
                    push!(ps_bs, p_bs)
                end

            end
            # Experimental
            vline!(p, [Kd_av]; label = "", lw = 6, color = :black, alpha = 0.3)
            vline!(p_bs, [Kd_av]; label = "", lw = 6, color = :black, alpha = 0.3)
            
            plot!(p; xlim = [m - margin, M + margin], size)
            plot!(p_bs; xlim = [m - margin, M + margin], size)
            push!(ps, p)
        end

        for k in [:xi, :D, :sGLC]
            p = plot(;title = Kd_ider, size)
            xticks =  (EXPS, string.(EXPS))
            vals = [Kd.val(k, exp) for exp in EXPS]
            p = bar!(p, EXPS, vals; title = k, label = "", xticks)
            push!(ps, p)
            push!(ps_bs, p)
        end

        pname = string(Kd_ider, "_marginals")
        mysavefig(ps, pname)

        pname = string(Kd_ider, "_marginals_vs_beta")
        mysavefig(ps_bs, pname; method2)
    end

end 

## -------------------------------------------------------------------
# marginals v2
let 
    objider = iJR.KAYSER_BIOMASS_IDER
    size = [300, 250]

    # Iders
    model_iders, Kd_iders = [objider], ["D"]
    for Kd_met in CONC_IDERS
        model_met = Kd_mets_map[Kd_met]
        model_exch = Kd_rxns_map[Kd_met]
        push!(model_iders, model_exch)
        push!(Kd_iders, string("u", Kd_met))
    end
    
    for (model_ider, Kd_ider) in zip(model_iders, Kd_iders)
        marg_params = (;xlabel = string(Kd_ider), yaxis = nothing, ylabel = "prob")

        epps = Plots.Plot[]
        exps = Plots.Plot[]
        for method in ALL_METHODS
            expp = plot(;title = string("Experimental"), marg_params...)
            epp = plot(;title = string(" MaxEnt: ", method), marg_params...)
            margin, m, M = -Inf, Inf, -Inf
            
            # EP
            for exp in EXPS
                Kd_av = Kd.val(Kd_ider, exp)
                color = exp_colors[exp]    

                datfile = INDEX[method, :DFILE, exp]
                dat = deserialize(datfile)
                model = dat[:model]
                objidx = ChU.rxnindex(model, objider)
                epouts = dat[:epouts]
                exp_beta = maximum(keys(epouts))
                epout = epouts[exp_beta]
                ep_av = ChU.av(model, epout, model_ider)
                ep_va = sqrt(ChU.va(model, epout, model_ider))
                        
                ChP.plot_marginal!(epp, model, [epout], model_ider; 
                    legend = false, color, alpha = 0.8, lw = 3)
                
                m = minimum([m, ep_av, Kd_av])
                M = maximum([M, ep_av, Kd_av])
                margin = maximum([margin, 3 * ep_va])

                # Experimental
                vline!(expp, [Kd_av]; label = "", lw = 3, color, alpha = 0.8)
                
            end
            
            map([expp, epp]) do p
                plot!(p; xlim = [m - margin, M + margin], size)
            end

            push!(epps, epp)
            push!(exps, expp)
        end

        extras = Plots.Plot[]
        for k in [:xi, :D, :sGLC]
            p = plot(;title = "Experimental", size, 
                xlabel = "rep", ylabel = string(k))
            xticks =  (EXPS, string.(EXPS))
            vals = [Kd.val(k, exp) for exp in EXPS]
            color = [exp_colors[exp] for exp in EXPS]
            p = bar!(p, EXPS, vals; label = "", xticks, color)
            push!(extras, p)
        end

        ps = Plots.Plot[exps; epps; extras]
        layout = (3, 3)
        pname = string(Kd_ider, "_marginals_v2")
        mysavefig(ps, pname; layout)

    end # for (model_ider, Kd_ider)

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