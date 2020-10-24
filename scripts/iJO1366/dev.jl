import DrWatson: quickactivate
    quickactivate(@__DIR__, "Chemostat_Kayser2005")

import Chemostat_Kayser2005: Chemostat
import Chemostat_Kayser2005.Chemostat.LP: MathProgBase
const Ch = Chemostat
const ChU = Chemostat.Utils
const ChSS = Chemostat.SteadyState
const ChLP = Chemostat.LP
const ChSU = Chemostat.SimulationUtils

import Chemostat_Kayser2005: KayserData, iJO1366
const iJO = iJO1366
const Kd = KayserData

## -------------------------------------------------------
const cache_dir = iJO.MODEL_CACHE_DATA_DIR
ChU.set_cache_dir(cache_dir)

## -------------------------------------------------------
sim_global_id = "MAXENT_FBA_EP_v1"
bundles = Dict()
for (Di, D) in Kd.val("D") |> enumerate
    ## SIMULATION PARAMS
    ξs = [Kd.val(:xi, Di)]
    βs = [0.0; ChU.logspace(6.0, 7.0, 10)]
    
    to_map = Iterators.product(ξs, [βs], Di)
    map(to_map) do (ξ, βs, Di)
        sim_hash = (sim_global_id, ξ, Di)
        res_id = (:RESULT, sim_hash)

        cached_dat = ChU.load_cache(res_id; verbose = false)
        isnothing(cached_dat) && return
        Di, ξ, βs, model, dat = cached_dat
        
        # boundling
        bundle = get!(bundles, Di, ChU.ChstatBundle())
        
        bundle[ξ, :net] = model |> ChU.compressed_model
        bundle[ξ, :fba] = dat[:fba]
        
        for (βi, β) in βs |> enumerate
            bundle[ξ, β, :ep] = dat[(:ep, βi)]
        end
        
    end
    return
end;

## -------------------------------------------------------
# SAVING
ChU.save_data(iJO.MAXENT_FBA_EB_BOUNDLES_FILE, bundles)
