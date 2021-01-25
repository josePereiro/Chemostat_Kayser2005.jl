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
    import JuMP.MathOptInterface
    using Serialization
    import UtilsJL
    const UJL = UtilsJL
    
    using Plots
    import GR
    GR.inline("png")
end

## -----------------------------------------------------------------------------------------------
LPDAT = ChU.load_data(iJR.LP_DAT_FILE)
const FBA_BOUNDED = :FBA_BOUNDEDs
const FBA_OPEN = :FBA_OPEN
const YIELD = :YIELD

## -----------------------------------------------------------------------------------------------
fileid = "3.1"
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
    FBA_OPEN => :red,
    FBA_BOUNDED => :orange,
    YIELD => :blue,
)

## -----------------------------------------------------------------------------------------------
# yield correlation
let
    p = plot(;title = "Yield correlation", xlabel = "exp", ylabel = "model")
    m, M = Inf, -Inf
    for exp in EXPS
        try
            model = LPDAT[YIELD, :model, exp]
            yout = LPDAT[YIELD, :yout, exp]
            model_yield = LPDAT[YIELD, :yield, exp]
            
            objidx = ChU.rxnindex(model, iJR.KAYSER_BIOMASS_IDER)
            exp_growth = Kd.val("D", exp)
            model_growth = ChU.av(model, yout, objidx)

            diff = abs(model_growth - exp_growth)/exp_growth
            diff > 0.05 && continue # unfeasible

            exp_yield = abs(Kd.val("D", exp) / Kd.uval("GLC", exp))
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
        p = plot(;title = "yield vs $(id)", xlabel = "exp $id", ylabel = "yield")
        for exp in EXPS
            try
                model = LPDAT[YIELD, :model, exp]
                yout = LPDAT[YIELD, :yout, exp]
                model_yield = LPDAT[YIELD, :yield, exp]
                # status, yflxs, model_yield, d, model = DAT[D]
                exglc_idx = ChU.rxnindex(model, "EX_glc_LPAREN_e_RPAREN__REV")
                biomass_idx = ChU.rxnindex(model, iJR.KAYSER_BIOMASS_IDER)
                
                Hd_val = Kd.val(id, exp)
                exp_yield = abs(Kd.val("D", exp) / Kd.uval("GLC", exp))
                scatter!(p, [Hd_val], [model_yield]; ms = 8,
                        color = :blue, alpha = 0.6, label = ""
                )
                scatter!(p, [Hd_val], [exp_yield]; ms = 8,
                        color = :red, alpha = 0.6, label = ""
                )
            catch err; @warn("Fail", err) end
        end
        push!(ps, p)
    end
    pname = "yield_vs_stuff"
    mysavefig(ps, pname)
end

## -----------------------------------------------------------------------------------------------
# let
#     for exp in Hd.EXPS
#         !haskey(DAT, exp) && continue
#         status, yflxs, yield, d, model = DAT[exp]
#         # @assert all(isapprox.(model.S * yflxs, model.b; atol = 1-3))
#         @show maximum(abs.(model.S * yflxs .- model.b))

#         break
#     end
# end

# ## -----------------------------------------------------------------------------------------------
# # correlations
# FLX_IDERS_MAP = Dict(
#     "GLC" => "EX_glc_LPAREN_e_RPAREN__REV",
#     "CO2" => "EX_co2_LPAREN_e_RPAREN_",
#     "O2" => "EX_o2_LPAREN_e_RPAREN__REV",
#     "AC" => "EX_ac_LPAREN_e_RPAREN_",
#     "NH4" => "EX_nh4_LPAREN_e_RPAREN__REV",
#     "D" => iJR.KAYSER_BIOMASS_IDER
# )

# COLORS = Dict(
#     "EX_glc_LPAREN_e_RPAREN__REV" => :blue,
#     "EX_co2_LPAREN_e_RPAREN_" => :yellow,
#     "EX_o2_LPAREN_e_RPAREN__REV" => :orange,
#     "EX_ac_LPAREN_e_RPAREN_" => :red,
#     "EX_nh4_LPAREN_e_RPAREN__REV" => :purple,
#     iJR.KAYSER_BIOMASS_IDER => :black,
# )

# let
#     yield_p = plot(title = "yield tot corrs"; xlabel = "exp flx", ylabel = "model flx")
#     fba_p = plot(title = "fba tot corrs"; xlabel = "exp flx", ylabel = "model flx")
#     margin, m, M = -Inf, Inf, -Inf
#     for (Kd_ider, model_ider) in FLX_IDERS_MAP
#         Hd_fun = Kd_ider == "D" ? Kd.val : Kd.uval
#         for (exp, D) in Kd.val(:D) |> enumerate
#             try
#                 !haskey(DAT, D) && continue
#                 status, yflxs, yield, d, model = DAT[D]
#                 model_idx = ChU.rxnindex(model, model_ider)
                
#                 color = COLORS[model_ider]
#                 Hd_flx = Hd_fun(Kd_ider, D)

#                 # yield
#                 ymax_flx = SENSE[model_ider] * yflxs[model_idx]
#                 scatter!(yield_p, [Hd_flx], [ymax_flx]; ms = 8,
#                     color, alpha = 0.6, label = ""
#                 )

#                 # fba
#                 fbaout = ChLP.fba(model, iJR.KAYSER_BIOMASS_IDER, iJR.COST_IDER)
#                 fba_flx = ChU.av(model, fbaout, model_ider)
#                 scatter!(fba_p, [Hd_flx], [fba_flx]; ms = 8,
#                     color, alpha = 0.6, label = ""
#                 )

#                 m = minimum([m, Hd_flx, ymax_flx, fba_flx])
#                 M = maximum([M, Hd_flx, ymax_flx, fba_flx])
#             catch err; @warn("Fail", err) end
#         end
#     end
#     margin = abs(M - m)*0.1
#     plot!(yield_p, [m - margin,M + margin], [m - margin,M + margin]; 
#         ls = :dash, color = :black, label = "")
#     plot!(fba_p, [m - margin,M + margin], [m - margin,M + margin]; 
#         ls = :dash, color = :black, label = "")
#     pname = "tot_corr"
#     mysavefig([yield_p, fba_p], pname)
# end

# ## -------------------------------------------------------------------
# # # leyends
# # # TODO fix this...
# # let
# #     for (title, colors) in [
# #             ("exp", exp_colors), 
# #             ("iders", ider_colors),
# #             ("method", method_colors)
# #         ]
# #     p = plot(; framestyle = :none)
# #         scolors = sort(collect(colors); by = (p) -> string(first(p)))
# #         for (id, color) in scolors
# #             scatter!(p, [0], [0];
# #                 thickness_scaling = 1,
# #                 color, ms = 8, label = string(id),
# #                 legendfontsize=10, 
# #                 # size = [300, 900],
# #                 # legend = :left
# #             )
# #         end
# #         mysavefig(p, "$(title)_color_legend")
# #     end
# # end