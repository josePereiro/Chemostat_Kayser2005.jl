# Data from Kayser2005 https://doi.org/10.1099/mic.0.27481-0.
module KayserData
    
    import Chemostat
    const Ch = Chemostat

    using ProjAssistant
    @gen_sub_proj

    include("data.jl")

    function __init__()
        _load_and_bundle()
        @create_proj_dirs
    end
end
