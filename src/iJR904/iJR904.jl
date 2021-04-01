module iJR904

    import BSON
    import ..Chemostat_Kayser2005
    const ChK = Chemostat_Kayser2005
    import ..BegData
    const Bd = BegData
    import Chemostat
    const ChU = Chemostat.Utils
    import ..KayserData
    const Kd = KayserData

    import UtilsJL
    const UJL = UtilsJL
    UJL.gen_sub_proj(@__MODULE__)


    include("const.jl")
    include("load_data.jl")
    include("beg_enz_cost.jl")
    include("kayser_biomass.jl")
    include("load_model.jl")

    function __init__()
        UJL.create_proj_dirs(@__MODULE__)
    end

end