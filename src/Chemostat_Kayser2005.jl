module Chemostat_Kayser2005
    
    import BSON
    import DrWatson
    import Chemostat
    const Ch = Chemostat
    const ChU = Chemostat.Utils
    const ChSS = Chemostat.SteadyState
    const ChLP = Chemostat.LP

    import UtilsJL
    const UJL = UtilsJL
    UJL.gen_top_proj(@__MODULE__)

    include("Utils/Utils.jl")
    include("BegData/BegData.jl")
    include("KayserData/KayserData.jl")
    include("iJR904/iJR904.jl")

    function __init__()
        UJL.create_proj_dirs(@__MODULE__)
    end

end # module
