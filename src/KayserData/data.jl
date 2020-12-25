## ------------------------------------------------------------------
# BUNDLE
const BUNDLE = Dict()
const MEDIUM = Dict()
const TABLE1 = Dict()
const TABLE2 = Dict()
function _load_and_bundle()

    empty!(MEDIUM)
    !isfile(KAYSER_CONV_MEDIUM_FILE) && return BUNDLE
    merge!(MEDIUM, load_data(KAYSER_CONV_MEDIUM_FILE; verbose=false))
    empty!(TABLE1)
    !isfile(KAYSER_CONV_TABLE1_FILE) && return BUNDLE
    merge!(TABLE1 ,load_data(KAYSER_CONV_TABLE1_FILE; verbose=false))
    !isfile(KAYSER_CONV_TABLE1_FILE) && return BUNDLE
    merge!(TABLE2 ,load_data(KAYSER_CONV_TABLE2_FILE; verbose=false))
    empty!(BUNDLE)

    for (id, dat) in TABLE1
        dict = get!(BUNDLE, id, Dict())
        dict["name"] = dat[1]
        dict["unit"] = dat[2]
        dict["vals"] = dat[3:end] .|> Float64
    end

    # for (id, dat) in TABLE2
    #     dict = get!(BUNDLE, id, Dict())
    #     dict["name"] = dat[1]
    #     dict["unit"] = dat[2]
    #     dict["vals"] = dat[3:end] .|> Float64
    # end

    let n = length(BUNDLE["sAC"]["vals"])
        for (id, dat) in MEDIUM
            dict = get!(BUNDLE, id, Dict())
            dict["name"] = dat["name"]
            dict["unit"] = dat["unit"]
            dict["vals"] = fill(float(dat["c"]), n)
        end
    end

    ## ------------------------------------------------------------------
    # xi
    dict = get!(BUNDLE, "xi", Dict())
    dict["name"] = "Cell specific dilution rate"
    dict["unit"] = "gDW/ L hr"
    dict["vals"] = BUNDLE["Xv"]["vals"] ./ BUNDLE["D"]["vals"];
    BUNDLE
end

## ------------------------------------------------------------------
# API
const msd_mets = ["AC", "GLC", "NH4"]
const iders_to_plot = ["AC", "GLC", "NH4", "D"]
val(dataid) = BUNDLE[string(dataid)]["vals"]
val(dataid, i) = val(dataid)[i]
cval(id) = val("c$id")
cval(id, i) = val("c$id", i)
sval(id) = val("s$id")
sval(id, i) = val("s$id", i)
name(dataid) = BUNDLE[string(dataid)]["name"]
unit(dataid) = BUNDLE[string(dataid)]["unit"]

