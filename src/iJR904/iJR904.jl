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

    using ProjAssistant
    @gen_sub_proj
    
    include("const.jl")
    include("load_data.jl")
    include("beg_enz_cost.jl")
    include("kayser_biomass.jl")
    include("load_model.jl")
    include("ME_MODES.jl")

    function __init__()
        @create_proj_dirs
    end

end