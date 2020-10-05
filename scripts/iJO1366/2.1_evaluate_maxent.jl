import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Kayser2005")

## ------------------------------------------------------------------

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
# GROWTH vs BETA
function growth_vs_beta()
    p = plot(title = nameof(iJO), 
        xlabel = "log beta", ylabel = "growth rate")
    for (Di, bundle) in bundles |> collect |> sort

        exp_growth = Kd.val("D", Di)
        plot!(p, log10.(bundle.βs), fill(exp_growth, length(bundle.βs)); 
            ls = :dash, color = colors[Di], lw = 3, label = "")
        
        exp_xi = Kd.val("xi", Di)
        fba_growth = ChU.av(bundle, exp_xi, :fba, iJO.OBJ_IDER)
        plot!(p, log10.(bundle.βs), fill(fba_growth, length(bundle.βs)); 
            ls = :dot, color = colors[Di], lw = 3, label = "")

        ep_growths = ChU.av(bundle, exp_xi, bundle.βs, :ep, iJO.OBJ_IDER)
        ep_stds = sqrt.(ChU.va(bundle, exp_xi, bundle.βs, :ep, iJO.OBJ_IDER))
        plot!(p, log10.(bundle.βs), ep_growths; 
            label = Di, lw = 3, color = colors[Di])
    end
    return p
end
growth_vs_beta()

## ------------------------------------------------------------------
# TOTAL STOI ERROR vs BETA
function total_stoi_err_vs_beta()
    p = plot(title = nameof(iJO), 
        xlabel = "log beta", ylabel = "stoi error [min/mean/max]")
    for (Di, bundle) in bundles |> collect |> sort
        exp_xi = Kd.val("xi", Di)
        model = bundle[exp_xi, :net]

        # collect errors
        M, N = size(model)
        errs = zeros(M, length(bundle.βs))
        for (βi, β) in bundle.βs |> enumerate
            epout = bundle[exp_xi, β, :ep]
            errs[:, βi] = ChU.norm1_stoi_err(model, epout)
        end

        # plots
        plot!(p, log10.(bundle.βs), maximum.(eachcol(errs));
            lw = 3, ls = :dash, c = colors[Di], label = "")
        plot!(p, log10.(bundle.βs), mean.(eachcol(errs));
            lw = 3, ls = :solid , c = colors[Di], label = string(Di))
        plot!(p, log10.(bundle.βs), minimum.(eachcol(errs));
            lw = 3, ls = :dot, c = colors[Di], label = "")
    end
    return p
end
total_stoi_err_vs_beta()

## ------------------------------------------------------------------
# METS STOI ERROR vs BETA
function mets_stoi_err_vs_beta(Di = 1)

    exch_met_map = iJO.load_exch_met_map()
    p = plot(title = string(nameof(iJO), " #", Di), 
        xlabel = "log beta", ylabel = "stoi error ")
    bundle = bundles[Di]
    exp_xi = Kd.val("xi", Di)
    model = bundle[exp_xi, :net]

    # collect errors
    err_dict = Dict()
    for (exch, _) in iJO.load_base_intake_info()
        met = exch_met_map[exch]
        for (βi, β) in bundle.βs |> enumerate
            epout = bundle[exp_xi, β, :ep]
            errs = get!(err_dict, met, [])
            push!(errs, ChU.norm1_stoi_err(model, epout, met))
        end
    end

    # plots
    serr_dict = sort(collect(err_dict); by = (p) -> sum(p[2]), rev = true)
    lcount = 5
    for (ider, errs) in serr_dict
        lb_ = lcount < 0 ? "" : ider
        plot!(p, log10.(bundle.βs), errs;
            lw = 3, ls = :dash, label = lb_)
        lcount -= 1
    end
    return p
end
mets_stoi_err_vs_beta()
