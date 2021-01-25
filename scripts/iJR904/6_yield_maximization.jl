import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Kayser2005")

@time begin
    import SparseArrays
    import Base.Threads: @threads, threadid, SpinLock

    #  ----------------------------------------------------------------------------
    # Run add https://github.com/josePereiro/Chemostat_Kayser2005.jl in the Julia Pkg REPL to install the
    # package, then you must activate the package enviroment (see README)
    import Chemostat_Kayser2005
    const ChK = Chemostat_Kayser2005

    const iJR = ChK.iJR904
    const Kd = ChK.KayserData # experimental data
    const Bd = ChK.BegData    # cost data

    #  ----------------------------------------------------------------------------
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

    import JuMP, GLPK
    using Serialization
    import UtilsJL
    const UJL = UtilsJL
    
    using Plots
    import GR
    GR.inline("png")
end

## -----------------------------------------------------------------------------------------------
function get_model(exp)
    BASE_MODELS = ChU.load_data(iJR.BASE_MODELS_FILE);
    model_dict = BASE_MODELS["fva_models"][exp]
    model = ChU.MetNet(;model_dict...) |> ChU.uncompressed_model
end
fileid = "6"
function mysavefig(p, pname; params...) 
    fname = UJL.mysavefig(p, string(fileid, "_", pname), iJR.MODEL_FIGURES_DIR; params...)
    @info "Plotting" fname
end
myminmax(a::Vector) = (minimum(a), maximum(a))
CONC_IDERS = String["GLC", "AcA", "FA"]
FLX_IDERS = [CONC_IDERS; "D"]

exp_colors = let
    colors = Plots.distinguishable_colors(length(Kd.EXPS))
    Dict(exp => color for (exp, color) in zip(Kd.EXPS, colors))
end

ider_colors = let
    colors = Plots.distinguishable_colors(length(FLX_IDERS))
    Dict(id => color for (id, color) in zip(FLX_IDERS, colors))
end

## -----------------------------------------------------------------------------------------------
function yLP(S, b, lb, ub, c, d; 
        ϵ = 1e-5, 
        sense = JuMP.MOI.MAX_SENSE,
        solver = GLPK.Optimizer
    )

    lp_model = JuMP.Model(solver)
    M, N = size(S)

    # Variables
    y = JuMP.@variable(lp_model, y[1:N])
    r = JuMP.@variable(lp_model, r)

    # Constraints
    JuMP.@constraint(lp_model, S * y - r .* b .== 0.0)
    JuMP.@constraint(lp_model, d' * y == 1.0)
    JuMP.@constraint(lp_model, y - r .* lb .>= 0.0)
    JuMP.@constraint(lp_model, y - r .* ub .<= 0.0)
    JuMP.@constraint(lp_model, r >= ϵ)

    # objective
    JuMP.@objective(lp_model, sense, c' * y)

    # optimize
    JuMP.optimize!(lp_model)

    yval = JuMP.value.(y)
    rval = JuMP.value.(r)
    v = yval ./ rval # sol
    y = (c' * v) / (d' * v) # yield

    status = JuMP.termination_status(lp_model)
    status, v, y
end
yLP(net, d; kwargs...) = yLP(net.S, net.b, net.lb, net.ub, net.c, d; kwargs...)

## -----------------------------------------------------------------------------------------------

# for method in [HOMO, EXPECTED, BOUNDED]
DAT = UJL.DictTree()
let
    # objider = iJR.KAYSER_BIOMASS_IDER
    for (exp, D) in Kd.val(:D) |> enumerate
        try 
            # prepare model
            model = get_model(exp)
            M, N = size(model)
            model = ChU.MetNet(model; c = zeros(N))
            ChU.check_dims(model)
            ChU.invert_bkwds!(model)
            exglc_idx = ChU.rxnindex(model, "EX_glc_LPAREN_e_RPAREN__REV")
            exac_idx = ChU.rxnindex(model, "EX_ac_LPAREN_e_RPAREN_")
            biomass_idx = ChU.rxnindex(model, iJR.KAYSER_BIOMASS_IDER)
            exp_growth = Kd.val(:D, exp)
            dgrowth = 0.0
            # ChU.lb!(model, iJR.KAYSER_BIOMASS_IDER, exp_growth * (1 - dgrowth))
            ChU.ub!(model, iJR.KAYSER_BIOMASS_IDER, exp_growth * (1 + dgrowth))

            # fba
            fbaout = ChLP.fba(model, iJR.KAYSER_BIOMASS_IDER)
            fba_growth = ChU.av(model, fbaout, iJR.KAYSER_BIOMASS_IDER)

            # yield max
            model.c[biomass_idx] = 1.0
            d = zeros(N); 
            d[exglc_idx] = 1.0
            # d[exac_idx] = 1.0
            status, yflxs, yield = yLP(model, d)
            status != JuMP.MOI.OPTIMAL && @warn status
            ymax_growth = yflxs[biomass_idx]
            DAT[D] = (;status, yflxs, yield, d, model)

            @info("Yield Maximization", 
                exp, status, yield,
                fba_growth, ymax_growth, exp_growth
            ); println()

        catch err
            @warn("Error", err, exp); println()
        end
    end
end


## -----------------------------------------------------------------------------------------------
# yield correlation
let
    p = plot(;title = "Yield correlation", xlabel = "exp", ylabel = "model")
    m, M = Inf, -Inf
    for (exp, D) in Kd.val(:D) |> enumerate
        exp in [1,2] && continue
        !haskey(DAT, D) && continue
        try
            status, yflxs, model_yield, d, model = DAT[D]
            exglc_idx = ChU.rxnindex(model, "EX_glc_LPAREN_e_RPAREN__REV")
            biomass_idx = ChU.rxnindex(model, iJR.KAYSER_BIOMASS_IDER)
            
            exp_yield = 1
            exp_yield = abs(D / Kd.uval("GLC", D))
            scatter!(p, [exp_yield], [model_yield]; ms = 8,
                    color = :blue, alpha = 0.6, label = ""
            )
            m = minimum([m, exp_yield, model_yield])
            M = maximum([M, exp_yield, model_yield])
        catch err; @warn("Fail", err) end
    end
    plot!(p, [m,M], [m,M]; ls = :dash, color = :black, label = "")
    pname = "yield_corr"
    mysavefig(p, pname)
end

## -----------------------------------------------------------------------------------------------
# yield vs stuff
let
    ps = Plots.Plot[]
    for id in [:D, :cGLC, :xi, :uGLC]
        p = plot(;title = "yield vs $(id)", xlabel = "exp $id", ylabel = "modeled yield")
        for (exp, D) in Kd.val(:D) |> enumerate
            !haskey(DAT, D) && continue
            try
                status, yflxs, model_yield, d, model = DAT[D]
                exglc_idx = ChU.rxnindex(model, "EX_glc_LPAREN_e_RPAREN__REV")
                biomass_idx = ChU.rxnindex(model, iJR.KAYSER_BIOMASS_IDER)
                
                Hd_val = Kd.val(id, D)
                scatter!(p, [Hd_val], [model_yield]; ms = 8,
                        color = :blue, alpha = 0.6, label = ""
                )
            catch err; @warn("Fail", err) end
        end
        push!(ps, p)
    end
    pname = "yield_vs_stuff"
    mysavefig(ps, pname)
end

# -----------------------------------------------------------------------------------------------
let
    for exp in Hd.EXPS
        !haskey(DAT, exp) && continue
        status, yflxs, yield, d, model = DAT[exp]
        # @assert all(isapprox.(model.S * yflxs, model.b; atol = 1-3))
        @show maximum(abs.(model.S * yflxs .- model.b))

        break
    end
end

## -----------------------------------------------------------------------------------------------
# correlations
FLX_IDERS_MAP = Dict(
    "GLC" => "EX_glc_LPAREN_e_RPAREN__REV",
    "CO2" => "EX_co2_LPAREN_e_RPAREN_",
    "O2" => "EX_o2_LPAREN_e_RPAREN__REV",
    "AC" => "EX_ac_LPAREN_e_RPAREN_",
    "NH4" => "EX_nh4_LPAREN_e_RPAREN__REV",
    "D" => iJR.KAYSER_BIOMASS_IDER
)

COLORS = Dict(
    "EX_glc_LPAREN_e_RPAREN__REV" => :blue,
    "EX_co2_LPAREN_e_RPAREN_" => :yellow,
    "EX_o2_LPAREN_e_RPAREN__REV" => :orange,
    "EX_ac_LPAREN_e_RPAREN_" => :red,
    "EX_nh4_LPAREN_e_RPAREN__REV" => :purple,
    iJR.KAYSER_BIOMASS_IDER => :black,
)

let
    yield_p = plot(title = "yield tot corrs"; xlabel = "exp flx", ylabel = "model flx")
    fba_p = plot(title = "fba tot corrs"; xlabel = "exp flx", ylabel = "model flx")
    margin, m, M = -Inf, Inf, -Inf
    for (Kd_ider, model_ider) in FLX_IDERS_MAP
        Hd_fun = Kd_ider == "D" ? Kd.val : Kd.uval
        for (exp, D) in Kd.val(:D) |> enumerate
            try
                !haskey(DAT, D) && continue
                status, yflxs, yield, d, model = DAT[D]
                model_idx = ChU.rxnindex(model, model_ider)
                
                color = COLORS[model_ider]
                Hd_flx = Hd_fun(Kd_ider, D)

                # yield
                ymax_flx = SENSE[model_ider] * yflxs[model_idx]
                scatter!(yield_p, [Hd_flx], [ymax_flx]; ms = 8,
                    color, alpha = 0.6, label = ""
                )

                # fba
                fbaout = ChLP.fba(model, iJR.KAYSER_BIOMASS_IDER, iJR.COST_IDER)
                fba_flx = ChU.av(model, fbaout, model_ider)
                scatter!(fba_p, [Hd_flx], [fba_flx]; ms = 8,
                    color, alpha = 0.6, label = ""
                )

                m = minimum([m, Hd_flx, ymax_flx, fba_flx])
                M = maximum([M, Hd_flx, ymax_flx, fba_flx])
            catch err; @warn("Fail", err) end
        end
    end
    margin = abs(M - m)*0.1
    plot!(yield_p, [m - margin,M + margin], [m - margin,M + margin]; 
        ls = :dash, color = :black, label = "")
    plot!(fba_p, [m - margin,M + margin], [m - margin,M + margin]; 
        ls = :dash, color = :black, label = "")
    pname = "tot_corr"
    mysavefig([yield_p, fba_p], pname)
end

## -------------------------------------------------------------------
# # leyends
# # TODO fix this...
# let
#     for (title, colors) in [
#             ("exp", exp_colors), 
#             ("iders", ider_colors),
#             ("method", method_colors)
#         ]
#     p = plot(; framestyle = :none)
#         scolors = sort(collect(colors); by = (p) -> string(first(p)))
#         for (id, color) in scolors
#             scatter!(p, [0], [0];
#                 thickness_scaling = 1,
#                 color, ms = 8, label = string(id),
#                 legendfontsize=10, 
#                 # size = [300, 900],
#                 # legend = :left
#             )
#         end
#         mysavefig(p, "$(title)_color_legend")
#     end
# end