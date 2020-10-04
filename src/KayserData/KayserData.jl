# Data from Kayser2005 https://doi.org/10.1099/mic.0.27481-0.
module KayserData

    import ..Chemostat_Kayser2005: 
        PROJ_ROOT, RAW_DATA_DIR, PROCESSED_DATA_DIR
    import ..Chemostat_Kayser2005.Chemostat.Utils: 
        load_data, save_data

    include("dirs_and_files.jl")
    include("data.jl")

    function __init__()
        _create_dirs()
    end
end
