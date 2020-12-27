## ------------------------------------------------------------------
# BUNDLE
const BUNDLE = Dict()
const MEDIUM = Dict()
const TABLE1 = Dict()
const TABLE2 = Dict()
function _load_and_bundle()

    empty!(BUNDLE)
    empty!(MEDIUM)
    empty!(TABLE1)
    empty!(TABLE2)

    !isfile(KAYSER_CONV_MEDIUM_FILE) && return BUNDLE
    !isfile(KAYSER_CONV_TABLE1_FILE) && return BUNDLE
    !isfile(KAYSER_CONV_TABLE2_FILE) && return BUNDLE
    
    merge!(MEDIUM, load_data(KAYSER_CONV_MEDIUM_FILE; verbose=false))
    merge!(TABLE1 ,load_data(KAYSER_CONV_TABLE1_FILE; verbose=false))
    merge!(TABLE2 ,load_data(KAYSER_CONV_TABLE2_FILE; verbose=false))

    for table in [TABLE1, TABLE2]
        for (id, dat) in table
            dict = get!(BUNDLE, id, Dict())
            dict["name"] = dat[1]
            dict["unit"] = dat[2]
            dict["vals"] = dat[3:end] .|> Float64
            dict["D"] = table["D"][3:end] .|> Float64
        end
    end
    BUNDLE["D"]["vals"] = sort(union(TABLE1["D"][3:end], TABLE2["D"][3:end]))
    BUNDLE["D"]["D"] = BUNDLE["D"]["vals"]

    let n = length(BUNDLE["D"]["vals"])
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
    dict["vals"] = BUNDLE["Xv"]["vals"] ./ BUNDLE["Xv"]["D"];
    BUNDLE
end

## ------------------------------------------------------------------
# API
const msd_mets = ["AC", "GLC", "NH4"]
const iders_to_plot = ["AC", "GLC", "NH4", "D"]
val(dataid) = BUNDLE[string(dataid)]["vals"]
val(dataid, exp::Int) = val(dataid)[exp]
function val(dataid, exp)
    dat = BUNDLE[dataid]
    i = findfirst(dat["D"] .== exp)
    val(dataid, i)
end
val(dataid, exp, dflt) = try; return val(dataid, exp); catch err; (@warn(err); dflt) end

cval(id) = val("c$id")
cval(id, exp) = val("c$id", exp)
cval(id, exp, dflt) = val("c$id", exp, dflt)

sval(id) = val("s$id")
sval(id, exp) = val("s$id", exp)
sval(id, exp, dflt) = val("s$id", exp, dflt)

uval(id) = val("u$id")
uval(id, exp) = val("u$id", exp)
uval(id, exp, dflt) = val("u$id", exp, dflt)

name(dataid) = BUNDLE[string(dataid)]["name"]
unit(dataid) = BUNDLE[string(dataid)]["unit"]

