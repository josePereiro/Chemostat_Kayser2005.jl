## ------------------------------------------------------------------
# BUNDLE
const bundle = Dict()
const medium = Dict()
const table1 = Dict()
function _load_and_bundle()

    empty!(medium)
    merge!(medium, load_data(KAYSER_CONV_MEDIUM_FILE; verbose=false))
    empty!(table1)
    merge!(table1 ,load_data(KAYSER_CONV_TABLE1_FILE; verbose=false))
    empty!(bundle)

    for (id, dat) in table1
        dict = get!(bundle, id, Dict())
        dict["name"] = dat[1]
        dict["unit"] = dat[2]
        dict["vals"] = dat[3:end] .|> Float64
    end

    let n = length(bundle["sAC"]["vals"])
        for (id, dat) in medium
            dict = get!(bundle, id, Dict())
            dict["name"] = dat["name"]
            dict["unit"] = dat["unit"]
            dict["vals"] = fill(float(dat["c"]), n)
        end
    end

    ## ------------------------------------------------------------------
    # xi
    dict = get!(bundle, "xi", Dict())
    dict["name"] = "Cell specific dilution rate"
    dict["unit"] = "gDW/ L hr"
    dict["vals"] = bundle["Xv"]["vals"] ./ bundle["D"]["vals"];
    bundle
end

## ------------------------------------------------------------------
# API
const msd_mets = ["AC", "GLC", "NH4"]
const iders_to_plot = ["AC", "GLC", "NH4", "D"]
val(dataid) = bundle[string(dataid)]["vals"]
val(dataid, i) = val(dataid)[i]
cval(id) = val("c$id")
cval(id, i) = val("c$id", i)
sval(id) = val("s$id")
sval(id, i) = val("s$id", i)
name(dataid) = bundle[string(dataid)]["name"]
unit(dataid) = bundle[string(dataid)]["unit"]

