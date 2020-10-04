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
# GROWTH evolution
p = plot(xlabel = "log beta", ylabel = "growth rate")
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
p

## ------------------------------------------------------------------
iJO.
