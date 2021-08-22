module Chemostat_Kayser2005
    
    import BSON
    import DrWatson
    import Chemostat
    const Ch = Chemostat
    const ChU = Chemostat.Utils
    const ChSS = Chemostat.SteadyState
    const ChLP = Chemostat.LP

    using ProjAssistant
    @gen_top_proj

    include("Utils/Utils.jl")
    include("BegData/BegData.jl")
    include("KayserData/KayserData.jl")
    include("iJR904/iJR904.jl")

end # module
