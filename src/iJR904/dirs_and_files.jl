# DIRS
const MODEL_RAW_DIR = joinpath(RAW_DATA_DIR, "iJR904")
const MODEL_PROCESSED_DATA_DIR = joinpath(PROCESSED_DATA_DIR, basename(MODEL_RAW_DIR))
const MODEL_FIGURES_DIR = joinpath(FIGURES_DATA_DIR, basename(MODEL_RAW_DIR))
const MODEL_CACHE_DIR = joinpath(MODEL_PROCESSED_DATA_DIR, "cache")

function _create_dirs()
    for dir in [MODEL_PROCESSED_DATA_DIR, 
                MODEL_FIGURES_DIR, MODEL_CACHE_DIR]
        try; mkdir(dir)
            catch err
        end
    end
end

# FILES
const MODEL_RAW_MAT_FILE = joinpath(MODEL_RAW_DIR, "iJR904.mat")
const BASE_MODELS_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "base_models.bson")
const BASE_BETA0_EPOUT_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "base_beta0_epout.bson")
const EXCH_MET_MAP_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "exch_met_map.bson")
const MAXENT_FBA_EB_RES_FILE = joinpath(MODEL_PROCESSED_DATA_DIR, "maxent_fba_ep_res.bson")