import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Kayser2005")

# ------------------------------------------------------------------
using MAT

import Chemostat
const ChLP = Chemostat.LP
const ChSS = Chemostat.SteadyState
const ChSU = Chemostat.SimulationUtils
const ChU = Chemostat.Utils

import Chemostat_Kayser2005: KayserData, iJO1366
const iJO = iJO1366
const Kd = KayserData

import ProgressMeter: ProgressThresh, Progress, next!, finish!, update!
using Plots

## ------------------------------------------------------------------
# Biomass medium sensivility
let
    model = ChU.load_data(iJO.BASE_MODEL_FILE; verbose = false)
    xi = Kd.val("xi") |> maximum
    intake_info = iJO.load_base_intake_info()
    results = Dict()
    factors = 0.0:0.1:1.0
    prog = Progress(length(factors) * length(intake_info))
    for (exch, dat) in intake_info
        c0 = dat["c"]
        res = get!(results, exch, [])
        for f in factors
            dat["c"] = f * c0
            ChSS.apply_bound!(model, xi, intake_info; emptyfirst = true)
            growth = try 
                    ChU.av(model, ChLP.fba(model, iJO.BIOMASS_IDER, iJO.COST_IDER), iJO.BIOMASS_IDER)
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

    p = plot(title = "Biomass medium sensivility", 
            xlabel = "fraction of initial conc", 
            ylabel = "growth"
        )
    th = 0.2
    sresults = sort(collect(results); by = (x) -> sum(x[2]))
    lcount = 4
    for (exch, res) in sresults
        lb_ = lcount > 0 ? exch : ""
        plot!(p, factors, res; label = lb_, lw = 3)
        lcount -= 1
    end
    p
end

## ------------------------------------------------------------------
# Checking fba_obj_val < exp_obj_val
let to_map = Kd.val("D") |> enumerate
    for (Di, D) in to_map

        model = ChU.load_data(iJO.BASE_MODEL_FILE; verbose = false)

        ## Open intakes except Glucose
        ex_glc_ider = "EX_glc__D_e"
        intake_info = iJO.load_base_intake_info()
        intake_info[ex_glc_ider]["c"] = Kd.val(:cGLC, Di)

        # Reajust biomass
        
        
        # impose constraint
        xi = Kd.val(:xi, Di)
        ChSS.apply_bound!(model, xi, intake_info; emptyfirst = true)

        fbaout = ChLP.fba(model, iJO.BIOMASS_IDER, iJO.COST_IDER);
        fba_obj_val = ChU.av(model, fbaout, iJO.BIOMASS_IDER)
        fba_obj_val = ChU.av(model, fbaout, iJO.BIOMASS_IDER)
        fba_ex_glc_val = ChU.av(model, fbaout, ex_glc_ider)
        fba_ex_glc_b = ChU.bounds(model, ex_glc_ider)
        exp_obj_val = Kd.val("D", Di)
        ChU.tagprintln_inmw("FBA SOLUTION", 
            "\nxi:                      ", xi,
            "\nobj_ider:                ", iJO.BIOMASS_IDER,
            "\nfba fba_ex_glc_val:      ", fba_ex_glc_val,
            "\nfba fba_ex_glc_b:        ", fba_ex_glc_b,
            "\nfba obj_val:             ", fba_obj_val,
            "\nexp obj_val:             ", exp_obj_val,
            "\ncost_ider:               ", iJO.COST_IDER,
            "\nfba cost_val:            ", ChU.av(model, fbaout, iJO.COST_IDER),
            "\n\n"
        )
        (fba_obj_val < exp_obj_val) && @warn "fba objval < exp objval" fba_obj_val exp_obj_val
    end
end

## ------------------------------------------------------------------
# Find cGLC that fit experimental growth
let
    Di = 1
    cGLC = Kd.val(:cGLC, Di)
    Kd_growth = Kd.val(:D, Di)
    xi = Kd.val(:xi, Di)
    model = ChU.load_data(iJO.BASE_MODEL_FILE; verbose = false)
    intake_info = iJO.load_base_intake_info()

    function work_fun(cGLC)
        ## Open intakes except Glucose
        intake_info["EX_glc__D_e"]["c"] = first(cGLC)
        
        # impose constraint
        ChSS.apply_bound!(model, xi, intake_info; emptyfirst = true)

        ## fba
        fbaout = ChLP.fba(model, iJO.BIOMASS_IDER, iJO.COST_IDER)
        fba_growth = ChU.av(model, fbaout, iJO.BIOMASS_IDER)
        return [fba_growth]
    end

    ChSU.grad_desc(work_fun; x0 = [cGLC], x1 = [cGLC * 0.9], th = 1e-5, 
        C = [cGLC * 0.1], target = [Kd_growth], maxiters = 1000)

end

## ------------------------------------------------------------------
# Testing scaled model
let
    model = ChU.load_data(iJO.BASE_MODEL_FILE; verbose = false)
    fbaout = ChLP.fba(model, iJO.BIOMASS_IDER, iJO.COST_IDER);
    ChU.tagprintln_inmw("FBA SOLUTION", 
        "\nobj_ider:         ", iJO.BIOMASS_IDER,
        "\nsize:             ", size(model),
        "\nfba obj_val:      ", ChU.av(model, fbaout, iJO.BIOMASS_IDER),
        "\ncost_ider:        ", iJO.COST_IDER,
        "\nfba cost_val:     ", ChU.av(model, fbaout, iJO.COST_IDER),
        "\n\n"
    )
    model = ChU.well_scaled_model(model, 100.0; verbose = false)
    fbaout = ChLP.fba(model, iJO.BIOMASS_IDER, iJO.COST_IDER);
    ChU.tagprintln_inmw("FBA SOLUTION", 
        "\nobj_ider:         ", iJO.BIOMASS_IDER,
        "\nsize:             ", size(model),
        "\nfba obj_val:      ", ChU.av(model, fbaout, iJO.BIOMASS_IDER),
        "\ncost_ider:        ", iJO.COST_IDER,
        "\nfba cost_val:     ", ChU.av(model, fbaout, iJO.COST_IDER),
        "\n\n"
    )
end