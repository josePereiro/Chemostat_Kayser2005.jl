using DrWatson 
quickactivate(@__DIR__, "Chemostat_Kayser2005")

using MAT

import Chemostat
import Chemostat.LP: fba, fva_preprocess
import Chemostat.SteadyState: apply_bound!
import Chemostat.Utils: MetNet, to_symbol_dict, isrev, split_revs, rxn_mets,
                        rxnindex, metindex, exchanges, expanded_model,
                        Rxn, Met, findempty,
                        av, va, nzabs_range, set_met!, set_rxn!,
                        isfixxed, ub, ub!, lb!, lb, FWD_SUFFIX, BKWD_SUFFIX

import Chemostat_Kayser2005: KayserData, BegData
import Chemostat_Kayser2005.iJR904: MODEL_RAW_MAT_FILE, Kd_mets_map, beg_enz_cost, 
                                    OBJ_IDER, ATPM_IDER, COST_IDER, beg_rxns_map, 
                                    ABS_MAX_BOUND, MAX_CONC, EXCH_MET_MAP_FILE
const Kd = KayserData

import UtilsJL: save_data, load_data, tagprintln_inmw

## ------------------------------------------------------------------
# LOAD RAW MODEL
src_file = MODEL_RAW_MAT_FILE
mat_model = MAT.matread(src_file)["model"]
model = MetNet(mat_model; reshape=true)
tagprintln_inmw("MAT MODEL LOADED", 
    "\nfile:             ", relpath(src_file), 
    "\nfile size:        ", filesize(src_file), " bytes", 
    "\nmodel size:       ", size(model),
    "\nnzabs_range:      ", nzabs_range(model.S),
)

## -------------------------------------------------------------------
# the maximal experimental growth rate in Kayser2005 is ~0.4 1/h
# The raw model present a growth rate bigger than that, so it is ok
# to use it directly as base model
fbaout = fba(model, OBJ_IDER);
tagprintln_inmw("FBA SOLUTION", 
    "\nobj_ider:         " , OBJ_IDER,
    "\nfba obj_val:      ", av(model, fbaout, OBJ_IDER),
    "\nmax exp obj_val:  ", maximum(Kd.val("D"))
)

## -------------------------------------------------------------------
# Set bounds
# The abs maximum bounds will be set to 100
tagprintln_inmw("CLAMP BOUNDS", 
    "\nabs max bound: ", ABS_MAX_BOUND
)
foreach(model.rxns) do ider
        isfixxed(model, ider) && return # fixxed reaction are untouched

        old_ub = ub(model, ider)
        new_ub = old_ub == 0.0 ? 0.0 : ABS_MAX_BOUND
        ub!(model, ider, new_ub)

        old_lb = lb(model, ider)
        new_lb = old_lb == 0.0 ? 0.0 : -ABS_MAX_BOUND
        lb!(model, ider, new_lb)
end

## -------------------------------------------------------------------
# CLOSING EXCHANGES
exchs = exchanges(model)
tagprintln_inmw("CLOSE EXCANGES", 
    "\nexchanges: ", exchs |> length
)
# Close, for now, all exchanges for avoiding it to be in revs
# The reversible reactions will be splited for modeling cost
# Exchanges have not associated cost, so, we do not split them
foreach(exchs) do idx
    ub!(model, idx, 0.0) # Closing all outtakes
    lb!(model, idx, 0.0) # Closing all intakes
end

## -------------------------------------------------------------------
# EXCHANGE METABOLITE MAP
exch_met_map = Dict()
for rxni in exchs
    rxn = model.rxns[rxni]
    metis = rxn_mets(model, rxni)
    mets = model.mets[metis]
    length(mets) != 1 && continue 
    met = mets |> first
    exch_met_map[rxn] = met
    exch_met_map[met] = rxn
end
save_data(EXCH_MET_MAP_FILE, exch_met_map)

## -------------------------------------------------------------------
# ENZYMATIC COST INFO
# The cost will be introduced as a reaction, we follow the same cost models as 
# Beg et al. (2007): https://doi.org/10.1073/pnas.0609845104.
# A new balance equations is then added:
#        Σ(rᵢ*costᵢ) + tot_cost = 0
#    Because the cost coefficients (costᵢ) < 0 (it resamble a reactant), the system must allocate 
#    the fluxes (rᵢ) so that Σ(rᵢ*costᵢ) = tot_cost, and tot_cost
#    are usually bounded [0.0, 1.0]
cost_info = Dict()
fwd_ider(rxn) = string(rxn, FWD_SUFFIX);
bkwd_ider(rxn) = string(rxn, BKWD_SUFFIX);
for rxn in model.rxns
    # The exchanges, the atpm and the biomass are synthetic reactions, so, 
    # they have should not have an associated enzimatic cost 
    any(startswith.(rxn, ["EX_", "DM_"])) && continue
    rxn == OBJ_IDER && continue
    rxn == ATPM_IDER && continue
        
    # Only the internal, non reversible reactions have an associated cost
    # We will split the rev reactions, so we save the cost for both versions (fwd, bkwd)
    if isrev(model, rxn)
        cost_info[fwd_ider(rxn)] = -beg_enz_cost(rxn)
        cost_info[bkwd_ider(rxn)] = -beg_enz_cost(rxn)
    else
        cost_info[rxn] = -beg_enz_cost(rxn)
    end
end

## -------------------------------------------------------------------
# SPLITING REVS
tagprintln_inmw("SPLITING REVS", 
    "\nfwd_suffix:      ", FWD_SUFFIX,
    "\nbkwd_suffix:     ", BKWD_SUFFIX,
)
model = split_revs(model;

    get_fwd_ider = fwd_ider,
    get_bkwd_ider = bkwd_ider,
);

## -------------------------------------------------------------------
# ADDING COST REACCION
cost_met_id = "cost"
cost_exch_id = COST_IDER
tagprintln_inmw("ADDING COST", 
    "\ncosts to add: ", cost_info |> length,
    "\nmin abs coe:  ", cost_info |> values .|> abs |> minimum,
    "\nmax abs coe:  ", cost_info |> values .|> abs |> maximum,
    "\ncost met id:  ", cost_met_id,
    "\ncost exch id: ", cost_exch_id
)

M, N = size(model)
cost_met = Met(cost_met_id, S = collect(values(cost_info)), rxns = collect(keys(cost_info)), b = 0.0)
model = expanded_model(model, M + 1, N + 1)
set_met!(model, findempty(model, :mets), cost_met)
cost_exch = Rxn(cost_exch_id, S = [1.0], mets = [cost_met_id], lb = -ABS_MAX_BOUND, ub = 0.0, c = 0.0)
set_rxn!(model, findempty(model, :rxns), cost_exch);

## ------------------------------------------------------------------
# BASE INTAKE INFO
# from fba
intake_info = Dict(
    "EX_glc_LPAREN_e_RPAREN_" => Dict("c" => MAX_CONC, "lb" => -ABS_MAX_BOUND),
    "EX_nh4_LPAREN_e_RPAREN_" => Dict("c" => MAX_CONC, "lb" => -ABS_MAX_BOUND),
    "EX_o2_LPAREN_e_RPAREN_" => Dict("c" => MAX_CONC, "lb" => -ABS_MAX_BOUND),
    "EX_pi_LPAREN_e_RPAREN_" => Dict("c" => MAX_CONC, "lb" => -ABS_MAX_BOUND),
    "EX_so4_LPAREN_e_RPAREN_" => Dict("c" => MAX_CONC, "lb" => -ABS_MAX_BOUND),
)

# from paper
kayser_medium = load_data(Kd.KAYSER_CONV_MEDIUM_FILE; verbose = false)
for (id, dat) in kayser_medium
    kmet = id[2:end] # drop the 'c'
    mmet = Kd_mets_map[kmet]
    exch = exch_met_map[mmet]
    intake_info[exch] = Dict("c" => dat["c"], "lb" => -ABS_MAX_BOUND)
end

## -------------------------------------------------------------------
# SET BASE EXCHANGE
tagprintln_inmw("SETTING EXCHANGES") 
# To control the intakes just the metabolites defined in the 
# base_intake_info (The minimum medium) will be opened.
# The base model will be constraint as in a cultivation with 
# experimental minimum xi
# see Cossios paper (see README)
ξ = minimum(Kd.val(:xi))
# println("Minimum medium: ", iJR.base_intake_info)
foreach(exchs) do idx
    ub!(model, idx, ABS_MAX_BOUND) # Opening all outakes
    lb!(model, idx, 0.0) # Closing all intakes
end

# see Cossios paper (see README) for details in the Chemostat bound constraint
apply_bound!(model, ξ, intake_info);

# tot_cost is the exchange that controls the bounds of the 
# enzimatic cost contraint, we bound it to [0, 1.0]
lb!(model, cost_exch_id, 0.0);
ub!(model, cost_exch_id, 1.0);

## -------------------------------------------------------------------
# FVA PREPROCESSING
model = fva_preprocess(model, 
#     eps = 1-9, # This avoid blocking totally any reaction
    verbose = true);

##
fbaout = fba(model, OBJ_IDER, COST_IDER);
tagprintln_inmw("FBA SOLUTION", 
    "\nobj_ider:         " , OBJ_IDER,
    "\nfba obj_val:      ", av(model, fbaout, OBJ_IDER),
    "\nmax exp obj_val:  ", maximum(Kd.val("D")),
    "\ncost_ider:        " , COST_IDER,
    "\nfba cost_val:     ", av(model, fbaout, COST_IDER),
    "\n"
)
Chemostat.Utils.summary(model, fbaout)