import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Kayser2005")

@time begin
    import SparseArrays
    import Base.Threads: @threads, threadid, SpinLock

    using Serialization

    import Chemostat_Kayser2005
    const ChK = Chemostat_Kayser2005

    const iJR = ChK.iJR904
    const Kd = ChK.KayserData # experimental data
    const Bd = ChK.BegData    # cost data

    import Chemostat
    import Chemostat.LP: MathProgBase
    const Ch = ChK.Chemostat
    const ChU = Ch.Utils
    const ChSS = Ch.SteadyState
    const ChLP = Ch.LP
    const ChEP = Ch.MaxEntEP
    const ChSU = Ch.SimulationUtils

    import UtilsJL
    const UJL = UtilsJL
    using Serialization
    using Base.Threads
    UJL.set_cache_dir(iJR.MODEL_CACHE_DIR)
end

## -------------------------------------------------------------------
# globals
const WLOCK = ReentrantLock()
const SIM_GLOBAL_ID = "iJR904_MAXENT_VARIANTS"
const DAT_FILE_PREFFIX =  "maxent_ep_dat_"

const INDEX = UJL.DictTree()
function dat_file(name; kwargs...)
    fname = UJL.mysavename(name, "jls"; kwargs...)
    joinpath(iJR.MODEL_PROCESSED_DATA_DIR, fname)
end

## -------------------------------------------------------------------
function base_model(exp)
    BASE_MODELS = ChU.load_data(iJR.BASE_MODELS_FILE; verbose = false);
    model_dict = BASE_MODELS["fva_models"][exp]
    ChU.MetNet(;model_dict...) |> ChU.uncompressed_model
end

# -------------------------------------------------------------------
# METHOD VARIANTS
const ME_Z_OPEN_G_OPEN          = :ME_Z_OPEN_G_OPEN           # Do not use extra constraints
const ME_MAX_POL                = :ME_MAX_POL    # Match ME and Dy biom average and constraint av_ug
const ME_Z_EXPECTED_G_BOUNDED   = :ME_Z_EXPECTED_G_BOUNDED    # Match ME and Dy biom average and constraint av_ug
const ME_Z_EXPECTED_G_MOVING    = :ME_Z_EXPECTED_G_MOVING     # 
const ME_Z_FIXXED_G_BOUNDED     = :ME_Z_FIXXED_G_BOUNDED      # Fix biom around observed

## -------------------------------------------------------------------
# ME_Z_EXPECTED_G_MOVING
include("2.0.1_ME_Z_EXPECTED_G_MOVING.jl")

## -------------------------------------------------------------------
# ME_Z_EXPECTED_G_BOUNDED
include("2.0.2_ME_Z_EXPECTED_G_BOUNDED.jl")

## -------------------------------------------------------------------
# ME_Z_FIXXED_G_BOUNDED
include("2.0.3_ME_Z_FIXXED_G_BOUNDED.jl")

## -------------------------------------------------------------------
# ME_Z_OPEN_G_OPEN
include("2.0.4_ME_Z_OPEN_G_OPEN.jl")

## -------------------------------------------------------------------
# ME_MAX_POL
include("2.0.5_ME_MAX_POL.jl")

## -------------------------------------------------------------------
# save index
ChU.save_data(iJR.MAXENT_VARIANTS_INDEX_FILE, INDEX; verbose = false)