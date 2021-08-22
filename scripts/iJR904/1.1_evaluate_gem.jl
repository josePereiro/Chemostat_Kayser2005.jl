using ProjAssistant
@quickactivate 

# ------------------------------------------------------------------
@time begin
    using MAT
    using SparseArrays

    import Chemostat
    const ChLP = Chemostat.LP
    const ChSS = Chemostat.SteadyState
    const ChU = Chemostat.Utils

    import SimTools
    const SimT = SimTools

    import Chemostat_Kayser2005: KayserData, iJR904
    const iJR = iJR904
    const Kd = KayserData

    import ProgressMeter: ProgressThresh, Progress, next!, finish!, update!
    using Plots
    import GR
    !isinteractive() && GR.inline("png")
end

## ------------------------------------------------------------------
# Biomass medium sensibility
let
    model = iJR.load_model("base_model")
    obj_ider = iJR.BIOMASS_IDER
    xi = Kd.val("xi") |> minimum
    intake_info = iJR.load_base_intake_info()
    results = Dict()
    factors = 0.0:0.01:1.0
    prog = Progress(length(factors) * length(intake_info))
    for (exch, dat) in intake_info
        c0 = dat["c"]
        res = get!(results, exch, [])
        for f in factors
            dat["c"] = f * c0
            ChSS.apply_bound!(model, xi, intake_info; emptyfirst = true)
            growth = try 
                    fbaout = ChLP.fba(model, obj_ider, iJR.COST_IDER)
                    ChU.av(model, fbaout, obj_ider)
                catch err; 0.0 end
            push!(res, growth)

            next!(prog; showvalues = [
                    (:exch, exch),
                    (:f, f),
                    (:c0, c0),
                    (:c, dat["c"]),
                    (:growth, growth)
                ]
            )
        end 
        dat["c"] = c0 
    end
    finish!(prog)

    p = plot(
            title = "Biomass medium sensivility", 
            xlabel = "fraction of initial conc", 
            ylabel = "growth"
        )
    sresults = sort(collect(results); by = (x) -> sum(x[2]))
    lcount = 4
    for (exch, res) in sresults
        lb_ = lcount > 0 ? exch : ""
        plot!(p, factors, res; label = lb_, lw = 3)
        lcount -= 1
    end
    sfig(iJR, p, 
        @fileid, "Biomass_medium_sensibility", ".png"
    )
end

## ------------------------------------------------------------------
# Checking fba_obj_val < exp_obj_val
let 
    to_map = Kd.val("D") |> enumerate
    for (Di, D) in to_map

        model = iJR.load_model("base_model")

        ## Open intakes except Glucose
        intake_info = iJR.load_base_intake_info()
        intake_info[iJR.EX_GLC_IDER]["c"] = Kd.val(:cGLC, Di)
        
        # impose constraint
        xi = Kd.val(:xi, Di)
        ChSS.apply_bound!(model, xi, intake_info; emptyfirst = true)

        fbaout = ChLP.fba(model, iJR.BIOMASS_IDER, iJR.COST_IDER);
        fba_obj_val = ChU.av(model, fbaout, iJR.BIOMASS_IDER)
        fba_obj_val = ChU.av(model, fbaout, iJR.BIOMASS_IDER)
        fba_ex_glc_val = ChU.av(model, fbaout, iJR.EX_GLC_IDER)
        fba_ex_glc_b = ChU.bounds(model, iJR.EX_GLC_IDER)
        exp_obj_val = Kd.val("D", Di)

        println("FBA SOLUTION", 
            "\nxi:                      ", xi,
            "\nobj_ider:                ", iJR.BIOMASS_IDER,
            "\nfba fba_ex_glc_val:      ", fba_ex_glc_val,
            "\nfba fba_ex_glc_b:        ", fba_ex_glc_b,
            "\nfba obj_val:             ", fba_obj_val,
            "\nexp obj_val:             ", exp_obj_val,
            "\ncost_ider:               ", iJR.COST_IDER,
            "\nfba cost_val:            ", ChU.av(model, fbaout, iJR.COST_IDER),
            "\n\n"
        )
        (fba_obj_val < exp_obj_val) && @warn "fba objval < exp objval" fba_obj_val exp_obj_val
    end
end
