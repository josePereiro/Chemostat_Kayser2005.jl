# Beg et al. (2007): https://doi.org/10.1073/pnas.0609845104.
module BegData

    using DataFrames
    import ..Chemostat_Kayser2005: PROJ_ROOT, RAW_DATA_DIR, PROCESSED_DATA_DIR
    import CSV

    include("dirs_and_files.jl")
    include("enz_data.jl")

    function __init__()
        load_enz_data()
    end
    
end