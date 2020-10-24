import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Kayser2005")

## ------------------------------------------------------------------
using MAT

import Chemostat
import Chemostat.LP: fba, fva_preprocess, fva
import Chemostat.SteadyState: apply_bound!
import Chemostat.Utils: MetNet, to_symbol_dict, isrev, split_revs, rxn_mets,
                        rxnindex, metindex, exchanges, expanded_model,
                        Rxn, Met, findempty, summary,
                        av, va, nzabs_range, set_met!, set_rxn!,
                        isfixxed, ub, ub!, lb!, lb, FWD_SUFFIX, BKWD_SUFFIX,
                        save_data, load_data, tagprintln_inmw, well_scaled_model

import Chemostat_Kayser2005: KayserData, BegData, iJO1366
const iJO = iJO1366
const Kd = KayserData

import ProgressMeter: ProgressThresh, Progress, next!, finish!, update!
using Plots

## ------------------------------------------------------------------
model = load_data(iJO.BASE_MODEL_FILE; verbose = false)
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
        apply_bound!(model, xi, intake_info; emptyfirst = true)
        growth = try 
                av(model, fba(model, iJO.BIOMASS_IDER, iJO.COST_IDER), iJO.BIOMASS_IDER)
            catch err 
                0.0
            end
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

## ------------------------------------------------------------------
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
    global lcount -= 1
end
p

## ------------------------------------------------------------------
to_map = Kd.val("D") |> enumerate
for (Di, D) in to_map

    model = load_data(iJO.BASE_MODEL_FILE; verbose = false)

    ## Open intakes except Glucose
    intake_info = iJO.load_base_intake_info()
    intake_info["EX_glc__D_e"]["c"] = Kd.val(:cGLC, Di) * 1.2
    
    # impose constraint
    xi = Kd.val(:xi, Di)
    apply_bound!(model, xi, intake_info; emptyfirst = true)

    ## fba
    fbaout = fba(model, iJO.BIOMASS_IDER, iJO.COST_IDER);
    tagprintln_inmw("FBA SOLUTION", 
        "\nxi:               ", xi,
        "\nobj_ider:         ", iJO.BIOMASS_IDER,
        "\nfba obj_val:      ", av(model, fbaout, iJO.BIOMASS_IDER),
        "\nexp obj_val:      ", Kd.val("D", Di),
        "\ncost_ider:        ", iJO.COST_IDER,
        "\nfba cost_val:     ", av(model, fbaout, iJO.COST_IDER),
        "\n\n"
    )
    # break;
end

## ------------------------------------------------------------------
# Find cGLC that fit experimental growth
Di = 1
cGLC = Kd.val(:cGLC, Di)
Kd_growth = Kd.val(:D, Di)
xi = Kd.val(:xi, Di)
th = 1e-5
last_diff = 0.0
sense = -1
step = 10.0
prog = ProgressThresh(th)
model = load_data(iJO.BASE_MODEL_FILE; verbose = false)
intake_info = iJO.load_base_intake_info()
for it in 1:50

    ## Open intakes except Glucose
    intake_info["EX_glc__D_e"]["c"] = cGLC
    
    # impose constraint
    apply_bound!(model, xi, intake_info; emptyfirst = true)

    ## fba
    fbaout = fba(model, iJO.BIOMASS_IDER, iJO.COST_IDER)
    fba_growth = av(model, fbaout, iJO.BIOMASS_IDER)
    
    diff = abs(fba_growth - Kd_growth)/abs(Kd_growth)
    diff < th && break # Converged
    diff >= last_diff && (global sense *= -1) # If not getting better, change sense
    global cGLC += sense * step * diff
    
    update!(prog, diff; showvalues = [
            (:it, it),
            (:sense, sense),
            (:fba_growth, fba_growth),
            (:Kd_growth, Kd_growth),
            (:diff, diff),
            (:last_diff, last_diff),
            (:cGLC, cGLC),
        ]
    )

    global last_diff = diff
end
finish!(prog)

## ------------------------------------------------------------------
model = load_data(iJO.BASE_MODEL_FILE; verbose = false)
fbaout = fba(model, iJO.BIOMASS_IDER, iJO.COST_IDER);
tagprintln_inmw("FBA SOLUTION", 
    "\nobj_ider:         ", iJO.BIOMASS_IDER,
    "\nsize:             ", size(model),
    "\nfba obj_val:      ", av(model, fbaout, iJO.BIOMASS_IDER),
    "\ncost_ider:        ", iJO.COST_IDER,
    "\nfba cost_val:     ", av(model, fbaout, iJO.COST_IDER),
    "\n\n"
)
model = well_scaled_model(model, 100.0; verbose = false)
fbaout = fba(model, iJO.BIOMASS_IDER, iJO.COST_IDER);
tagprintln_inmw("FBA SOLUTION", 
    "\nobj_ider:         ", iJO.BIOMASS_IDER,
    "\nsize:             ", size(model),
    "\nfba obj_val:      ", av(model, fbaout, iJO.BIOMASS_IDER),
    "\ncost_ider:        ", iJO.COST_IDER,
    "\nfba cost_val:     ", av(model, fbaout, iJO.COST_IDER),
    "\n\n"
)

## ------------------------------------------------------------------
@show cGLC
@show Kd.val(:cGLC, Di) * 1.2
@show Kd.val(:cGLC, Di)
##
Kd.val(:cGLC, Di)/cGLC
##
model = load_data(iJO.BASE_MODEL_FILE; verbose = false)
lb!(model, iJO.COST_IDER, 0.0);
ub!(model, iJO.COST_IDER, 1.0);
# model = load_data(iJO.BASE_MODEL_FILE; verbose = false);
save_data(iJO.BASE_MODEL_FILE, model)
# fbaout = fba(model, iJO.BIOMASS_IDER, iJO.COST_IDER);
# tagprintln_inmw("FBA SOLUTION", 
#     "\nobj_ider:         ", iJO.BIOMASS_IDER,
#     "\nfba obj_val:      ", av(model, fbaout, iJO.BIOMASS_IDER),
#     "\ncost_ider:        ", iJO.COST_IDER,
#     "\nfba cost_val:     ", av(model, fbaout, iJO.COST_IDER),
#     "\n\n"
# )
