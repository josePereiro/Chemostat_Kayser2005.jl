import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Kayser2005")

## ------------------------------------------------------------------

import SparseArrays
import Chemostat_Kayser2005: Chemostat
import Chemostat_Kayser2005.Chemostat.LP: MathProgBase
const Ch = Chemostat
const ChU = Chemostat.Utils
const ChSS = Chemostat.SteadyState
const ChLP = Chemostat.LP
import Chemostat_Kayser2005: KayserData, iJO1366
const iJO = iJO1366
const Kd = KayserData
using Distributions

import ProgressMeter: ProgressThresh, Progress, next!, finish!, update!
import Plots: plot, plot!, scatter, scatter!

## ------------------------------------------------------------------
# LOAD DATA
bundles = ChU.load_data(iJO.MAXENT_FBA_EB_BOUNDLES_FILE)

## ------------------------------------------------------------------
color_pool = [:orange, :blue, :red, :black, :violet, 
    :gray, :green, :brown, :magenta];
colors = Dict(Di => rand(color_pool) for Di in keys(bundles))

## ------------------------------------------------------------------
closest_βs = Dict()
for (Di, bundle) in bundles
    exp_ξ = Kd.val("xi", Di)
    exp_μ = Kd.val("D", Di)
    closest_βs[Di] = ChU.find_closest_beta(bundle, exp_ξ, exp_μ, iJO.OBJ_IDER)
end
closest_βs

## -------------------------------------------------------------------
# CLOSEST BETA VS BETA
# Just for checking that the experimental objective is inside the beta intervals
# and evaluate the 'experimental' beta approximation
function closest_beta()
    p = plot(title = nameof(iJO), 
        xlabel = "log beta", ylabel = "growth rate"
    )
    for (Di, bundle) in bundles
        
        # model
        exp_ξ = Kd.val("xi", Di)
        μs = ChU.av(bundle, exp_ξ, bundle.βs, :ep, iJO.OBJ_IDER)
        plot!(p, log10.(bundle.βs), μs, label = "", color = colors[Di], lw = 3)

        # exp
        scatter!(p, [log10(closest_βs[Di])], [Kd.val("D", Di)], label = "", 
            color = colors[Di], ms = 9)
    end
    return p
end
closest_beta()

## -------------------------------------------------------------------
# GROWTH CORRELATION
function growth_correlation()
    p = plot(title = nameof(iJO), 
        xlabel = "experimental growth rate", ylabel = "modeled growth rate"
    )

    fun(x) = x#log10(x)
    for (Di, bundle) in bundles
        
        exp_ξ = Kd.val("xi", Di)
        exp_β = closest_βs[Di]
        model_μ = ChU.av(bundle, exp_ξ, exp_β, :ep, iJO.OBJ_IDER)
        exp_μ = Kd.val("D", Di)
        scatter!(p, [fun(exp_μ)], [fun(model_μ)], label = "", color = colors[Di])

    end
    return p
end
growth_correlation()

## ------------------------------------------------------------------
# TEST DATA
model = bundles[1][1, :net];
out = bundles[1][1, :fba];
intake_info = iJO.load_base_intake_info();

## ------------------------------------------------------------------
# CONC CORRELATION
function conc_correlations()
    conc_dict = Dict()
    mets_map = iJO.load_mets_map()
    exch_met_map = iJO.load_exch_met_map()

    fbap = plot(title = "FBA CORRELATION")
    epp = plot(title = "EP CORRELATION")

    for (Di, bundle) in bundles
        
        dat = get!(conc_dict, Di, Dict())
        exp_ξ =  Kd.val("xi", Di)
        exp_β = closest_βs[Di]

        for kmet in Kd.msd_mets
            mrxn = exch_met_map[mets_map[kmet]]

            # Exp
            exp_ss = get!(dat, "exp", [])
            s = Kd.sval(kmet, Di)
            push!(exp_ss, s)

            # fba
            fba_ss = get!(dat, "fba", [])
            (s, serr) = ChU.medium_conc(bundle, exp_ξ, :fba, mrxn) 
            push!(fba_ss, s)
            

            # ep
            ep_ss = get!(dat, "ep", [])
            ep_serrs = get!(dat, "eperr", [])
            (s, serr) = ChU.medium_conc(bundle, exp_ξ, exp_β, :ep, mrxn) 
            push!(ep_ss, s)
            push!(ep_serrs, serr)
        end

        fun(x) = log10(x + 13-4)
        scatter!(fbap, fun.(dat["exp"]), fun.(dat["fba"]), label = "", color = colors[Di])
        scatter!(epp, fun.(dat["exp"]), fun.(dat["ep"]), 
            # yerr = fun.(dat["eperr"]); 
            label = "", color = colors[Di])

    end

    # conc_dict
    plot([fbap, epp]..., layout = 2)
end
conc_correlations()


## ------------------------------------------------------------------
# # TOTAL CORRELATION
# function total_correlation()
#     ider_map = Kd.load_rath_met_exch_map()
#     ticks = round.(collect(range(-0.4, 0.4; length = 4)), digits = 2)
#     lims = [minimum(ticks), maximum(ticks)]
#     fbap = plot(
#         ylim = lims, 
#         xlim = lims, 
#         xticks = ticks,
#         yticks = ticks,
#         xtickfontsize = 9,
#         ytickfontsize = 9,
#         title = "FBA correlation", 
#         xlabel = "exp", ylabel = "model"
#     )
#     plot!(fbap, x -> x, lims[1], lims[2]; label = "", color = :black)
#     # ticks = round.(collect(range(-0.4, 0.4; length = 4)), digits = 2)
#     # lims = [minimum(ticks), maximum(ticks)]
#     epp = plot(
#         ylim = lims, 
#         xlim = lims, 
#         xticks = ticks,
#         yticks = ticks,
#         xtickfontsize = 9,
#         ytickfontsize = 9,
#         title = "EP correlation", 
#         xlabel = "exp", ylabel = "model"
#     )
#     plot!(epp, x -> x, lims[1], lims[2]; label = "", color = :black)

#     xs = []
#     xerrs = []
#     fbays = []
#     epys = []
#     epyerrs = []
#     for (Di, bundle) in bundles
#         exp_ξ =  Kd.val("xi", Di)
#         exp_μ =  Kd.val("D", Di)
#         exp_β = closest_βs[Di]


#         # ChU.medium_conc(bundles[1], 1, 1, :ep, "EX_cit_e"; 
#         for kider in Kd.iders_to_plot
#             sense = kider == "D" ? 1 : -1 # TODO: package this
#             mider = ider_map[kider]
#             exp_av = sense * Kd.sval(kider, Di)
#             push!(xs, exp_av)
#             push!(xerrs, exp_err)

#             # FBA
#             mod_av = ChU.av(bundle, exp_ξ, :fba, mider)
#             push!(fbays, mod_av)
            
#             # EP
#             mod_av = ChU.av(bundle, exp_ξ, exp_β,:ep, mider)
#             push!(epys, mod_av)
#             mod_err = sqrt(ChU.va(bundle, exp_ξ, exp_β,:ep, mider))
#             push!(epyerrs, mod_err)
#         end
#     end

#     # FBA
#     scatter!(fbap, xs, fbays; xerr = xerrs, label = "", color = :black)

#     # EP
#     scatter!(epp, xs, epys; xerr = epys, yerr = epyerrs, label = "", color = :black)

#     return plot([fbap, epp]..., layout = 2, size = [800, 400],)
# end
# total_correlation()