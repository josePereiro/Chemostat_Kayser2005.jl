module iJR904

    import BSON
    import ..Chemostat_Kayser2005
    import ..Chemostat_Kayser2005: PROJ_ROOT, DATA_DIR, FIGURES_DATA_DIR, RAW_DATA_DIR, PROCESSED_DATA_DIR
    import ..BegData
    const Bd = BegData
    const ChU = Chemostat_Kayser2005.Chemostat.Utils
    import ..KayserData
    const Kd = KayserData

    include("const.jl")
    include("dirs_and_files.jl")
    include("load_data.jl")
    include("beg_enz_cost.jl")
    include("kayser_biomass.jl")

    function __init__()
        _create_dirs()
    end

end