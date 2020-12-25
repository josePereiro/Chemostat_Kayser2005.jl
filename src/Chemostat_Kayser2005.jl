module Chemostat_Kayser2005
    
    import BSON
    import DrWatson
    import Chemostat
    const Ch = Chemostat
    const ChU = Chemostat.Utils
    const ChSS = Chemostat.SteadyState
    const ChLP = Chemostat.LP

    include("Utils/Utils.jl")
    include("BegData/BegData.jl")
    include("KayserData/KayserData.jl")
    include("iJR904/iJR904.jl")
    include("iJO1366/iJO1366.jl")

end # module
