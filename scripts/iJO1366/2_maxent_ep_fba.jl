import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Kayser2005")

## ------------------------------------------------------------------
import Chemostat
import Chemostat.LP: fba, fva_preprocess
import Chemostat.SteadyState: apply_bound!
import Chemostat.SimulationUtils: cached_simulation
import Chemostat.Utils: MetNet, rxnindex, metindex, summary,
                        save_data, load_data, tagprintln_inmw
import Chemostat_Kayser2005: KayserData, BegData, iJO1366
const iJO = iJO1366
const Kd = KayserData

## ------------------------------------------------------------------
function prepare_model(xi)
    model = load_data(iJO.BASE_MODEL_FILE; verbose = false)
    intake_info = iJO.load_base_intake_info()
    apply_bound!(model, xi, intake_info)
    return model
end

##
prepare_model(1)

## ------------------------------------------------------------------
tagprintln_inmw("Bla")