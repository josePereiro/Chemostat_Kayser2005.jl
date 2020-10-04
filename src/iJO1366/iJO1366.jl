# Download link http://bigg.ucsd.edu/models/iJO1366
module iJO1366

    import BSON
    import ..Chemostat_Kayser2005: 
        PROJ_ROOT, DATA_DIR, FIGURES_DATA_DIR, 
        RAW_DATA_DIR, PROCESSED_DATA_DIR
    import ..BegData
    const Bd = BegData
    import ..KayserData
    const Kd = KayserData
    import ..Chemostat_Kayser2005.Chemostat.Utils: 
        load_data, save_data

    include("const.jl")
    include("dirs_and_files.jl")
    include("load_data.jl")
    include("beg_enz_cost.jl")

    function __init__()
        _create_dirs()
    end

end