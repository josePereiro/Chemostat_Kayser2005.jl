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
fileid = "2.1"
function mysavefig(p, pname; params...)
    pname = string(fileid, "_", pname)
    fname = UJL.mysavefig(p, pname, iJR.MODEL_FIGURES_DIR; params...)
    @info("Plotting", fname)
end
myminmax(a::Vector) = (minimum(a), maximum(a))

Kd_rxns_map = iJR.load_Kd_rxns_map()
Kd_mets_map = iJR.load_Kd_mets_map()

## -------------------------------------------------------------------
# Collect
include("4.0.1_collect_DAT.jl")

# ## -------------------------------------------------------------------
# EXPS = DAT[:EXPS]

# exp_colors = let
#     colors = Plots.distinguishable_colors(length(EXPS))
#     Dict(exp => color for (exp, color) in zip(EXPS, colors))
# end

# ider_colors = Dict(
#     "GLC" => :red, "CO2" => :yellow,
#     "O2" => :blue, "AC" => :orange, 
#     "NH4" => :green, "D" => :black,
# )

# method_colors = Dict(
#     ME_Z_OPEN_G_OPEN => :red,
#     ME_MAX_POL => :blue,
#     ME_Z_EXPECTED_G_BOUNDED => :orange,
#     ME_Z_EXPECTED_G_MOVING => :purple,
#     ME_Z_FIXXED_G_BOUNDED => :green,
# )


# ## -------------------------------------------------------------------
# # PLOTS
# # -------------------------------------------------------------------

# ## -------------------------------------------------------------------
# # polytope box volume
# include("4.0.2_pol_vox_volume.jl")

# ## -------------------------------------------------------------------
# # Var study
# include("4.0.3_var_study.jl")

# ## -------------------------------------------------------------------
# # MSE study
# include("4.0.4_MSE_study.jl")

# ## -------------------------------------------------------------------
# # proj 2D
# include("4.0.5_proj2d_study.jl")

# ## -------------------------------------------------------------------
# # Beta tendencies
# include("4.0.6_Beta_tendencies.jl")

# ## -------------------------------------------------------------------
# # flx correlations
# include("4.0.7_flx_correlations.jl")

# ## -------------------------------------------------------------------
# # inbounds study
# include("4.0.7_flx_correlations.jl")

# ## -------------------------------------------------------------------
# # marginals
# include("4.0.9_marginals.jl")

# ## -------------------------------------------------------------------
# # legends
# include("4.0.10_leyends.jl")