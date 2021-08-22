# Beg et al. (2007): https://doi.org/10.1073/pnas.0609845104.
module BegData

    using DataFrames
    import CSV

    using ProjAssistant
    @gen_sub_proj
    
    # The average crowding coefficient a was fitted to obtain the minimum square 
    # deviation between the measured and model predicted growth rates, 
    # resulting in a = 0.0040 ± 0.0005 h•g/mmol, in which g is grams dry weight. 
    # However, the maximum growth rates on glucose and glycerol are more consistent 
    # with a = 0.0031 ± 0.0001 h•g/mmol and a = 0.0053 ± 0.0001 h•g/mmol, respectively.
    const AVE_COST  = 0.0031 # I will use this to set the unknown costs


    # enzymatic costs from Beg et al. (2007): https://doi.org/10.1073/pnas.0609845104.
    function load_enz_data()
        datfile = rawdir(BegData, "beg2007___enzymatic_data.tsv")
        return CSV.read(datfile, DataFrame)
    end

    function __init__()
        @create_proj_dirs
    end
    
end