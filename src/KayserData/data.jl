## ------------------------------------------------------------------
# HAYSER_CONV_TABLE1_FILE, HAYSER_CONV_MEDIUM_FILE
medium = load_data(HAYSER_CONV_MEDIUM_FILE; verbose=false)
table1 = load_data(HAYSER_CONV_TABLE1_FILE; verbose=false)

## ------------------------------------------------------------------
# BUNDLE
bundle = Dict()
for (id, dat) in table1
    dict = get!(bundle, id, Dict())
    n = length(dat)
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
# API
msd_mets = ["AC", "GLC", "NH4"]
val(dataid) = bundle[string(dataid)]["vals"]
val(dataid, i) = val(dataid)[i]
name(dataid) = bundle[string(dataid)]["name"]
unit(dataid) = bundle[string(dataid)]["unit"]

## ------------------------------------------------------------------
# xi
bundle["xi"] = Dict()
bundle["xi"]["name"] = "Cell specific dilution rate"
bundle["xi"]["unit"] = "gDW/ L hr"
bundle["xi"]["vals"] = val("Xv") ./ val("D");
