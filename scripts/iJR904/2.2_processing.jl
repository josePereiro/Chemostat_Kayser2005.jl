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
include("2.2.1_collect_DAT.jl")

## -------------------------------------------------------------------
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


## -------------------------------------------------------------------
# PLOTS
# -------------------------------------------------------------------

## -------------------------------------------------------------------
# polytope box volume
include("2.2.2_pol_vox_volume.jl")

## -------------------------------------------------------------------
# Var study
include("2.2.3_var_study.jl")

## -------------------------------------------------------------------
# MSE study
include("2.2.4_MSE_study.jl")

## -------------------------------------------------------------------
# proj 2D
include("2.2.5_proj2d_study.jl")

## -------------------------------------------------------------------
# Beta tendencies
include("2.2.6_Beta_tendencies.jl")

## -------------------------------------------------------------------
# flx correlations
include("2.2.7_flx_correlations.jl")

## -------------------------------------------------------------------
# inbounds study
include("2.2.7_flx_correlations.jl")

## -------------------------------------------------------------------
# marginals
include("2.2.9_marginals.jl")

## -------------------------------------------------------------------
# legends
include("2.2.10_leyends.jl")