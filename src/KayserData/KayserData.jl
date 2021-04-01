# Data from Kayser2005 https://doi.org/10.1099/mic.0.27481-0.
module KayserData

    import ..Chemostat_Kayser2005
    const ChK = Chemostat_Kayser2005
    import Chemostat
    const Ch = Chemostat

    import UtilsJL
    const UJL = UtilsJL
    UJL.gen_sub_proj(@__MODULE__)

    include("data.jl")

    function __init__()
        _load_and_bundle()
        UJL.create_proj_dirs(@__MODULE__)
    end
end
