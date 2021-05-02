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

    CONV_MEDIUM_FILE = procdir("medium_conv.bson")
    CONV_TABLE1_FILE = procdir("table1_conv.bson")
    CONV_TABLE2_FILE = procdir("table2_conv.bson")

    !isfile(CONV_MEDIUM_FILE) && return BUNDLE
    !isfile(CONV_TABLE1_FILE) && return BUNDLE
    !isfile(CONV_TABLE2_FILE) && return BUNDLE
    
    merge!(MEDIUM, UJL.load_data(CONV_MEDIUM_FILE; verbose=false))
    merge!(TABLE1, UJL.load_data(CONV_TABLE1_FILE; verbose=false))
    merge!(TABLE2, UJL.load_data(CONV_TABLE2_FILE; verbose=false))

    # collect all in a single bundle
    for table in [TABLE1, TABLE2]
        for (id, dat) in table
            dict = get!(BUNDLE, id, Dict())
            dict["name"] = dat[1]
            dict["unit"] = dat[2]
            dict["vals"] = dat[3:end] .|> Float64
        end
    end
    BUNDLE["D"]["vals"] = TABLE1["D"][3:end]

    # the gases exchanges are non reported 
    # for the experiments [14, 15]...
    # I computed from an extrapolation with the CTR and OTR.
    for (uid, tid) in [("uCO2", "CTR"), ("uO2", "OTR")]
        u = BUNDLE[uid]["vals"] # mmol/ g hr
        tr = BUNDLE[tid]["vals"] # g/ L hr

        # reported indexes
        rep_idxs = 1:min(length(u), length(tr))

        # u = m * tr
        m = sum(u[rep_idxs] ./ tr[rep_idxs])/ length(rep_idxs)
        infer_u = m .* tr
        infer_u[rep_idxs] .= u[rep_idxs]
        @assert length(infer_u) == length(tr)

        BUNDLE[uid]["vals"] = infer_u
    end


    let n = length(BUNDLE["D"]["vals"])
        for (id, dat) in MEDIUM
            dict = get!(BUNDLE, id, Dict())
            dict["name"] = dat["name"]
            dict["unit"] = dat["unit"]
            dict["vals"] = fill(float(dat["c"]), n)
            dict["D"] = BUNDLE["D"]["vals"]
        end
    end

    ## ------------------------------------------------------------------
    # xi
    dict = get!(BUNDLE, "xi", Dict())
    dict["name"] = "Cell specific dilution rate"
    dict["unit"] = "gDW/ L hr"
    dict["vals"] = BUNDLE["Xv"]["vals"] ./ BUNDLE["D"]["vals"];

    ## ------------------------------------------------------------------
    # X (common name)
    BUNDLE["X"] = BUNDLE["Xv"]

    # exchage for exp 14 and 15 wans't directly reported
    # I compute it using 0 = uX + (c - s)D
    for met in ["AC", "GLC", "NH4"]
        for exp in [14, 15]
            
            # cval(met, exp, 0.0)
            cdat = get(BUNDLE, string("c", met), Dict("vals" => Dict(exp => 0.0)))
            c = cdat["vals"][exp]
            # sval(met, exp, 0.0)
            s = BUNDLE[string("s", met)]["vals"][exp] 
            # val(:D, exp)
            D = BUNDLE["D"]["vals"][exp]
            # val(:Xv, exp)
            Xv = BUNDLE["Xv"]["vals"][exp]
            u = -(c - s) * D / Xv

            udat = BUNDLE[string("u", met)]
            push!(udat["vals"], u)
        end
    end

    BUNDLE
end

## ------------------------------------------------------------------
# API
# We will drop exp 14 and 15 because 
# poor carbon recovery    
const EXPS = 1:13
const msd_mets = ["AC", "GLC", "NH4"]
const iders_to_plot = ["AC", "GLC", "NH4", "D"]
val(dataid) = BUNDLE[string(dataid)]["vals"][EXPS]
val(dataid, exp::Int) = val(dataid)[exp]
val(dataid, exp, dflt) = try; return val(dataid, exp); catch err; dflt end

cval(id, args...) = val("c$id", args...)
sval(id, args...) = val("s$id", args...)
uval(id, args...) = val("u$id", args...)

name(dataid) = BUNDLE[string(dataid)]["name"]
unit(dataid) = BUNDLE[string(dataid)]["unit"]

ciD_X(id) = [cval(id, exp, 0.0) * val(:D, exp) / val(:Xv, exp) for exp in EXPS]
ciD_X(id, exp) = cval(id, exp, 0.0) * val(:D, exp) / val(:Xv, exp)
