import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Kayser2005")

import SparseArrays
import Base.Threads: @threads, threadid, SpinLock

# -------------------------------------------------------------------
using Serialization

import Chemostat_Kayser2005
import Chemostat_Kayser2005.Chemostat
import Chemostat_Kayser2005.Chemostat.LP: MathProgBase
const ChK = Chemostat_Kayser2005

const iJR = ChK.iJR904
const Kd = ChK.KayserData # experimental data
const Bd = ChK.BegData    # cost data

const Ch = ChK.Chemostat
const ChU = Ch.Utils
const ChSS = Ch.SteadyState
const ChLP = Ch.LP
const ChEP = Ch.MaxEntEP
const ChSU = Ch.SimulationUtils

## -------------------------------------------------------------------
let
    # setup
    cache_dir = iJR.MODEL_CACHE_DIR

    # orig model
    exp = 4
    model = ChU.load_data(iJR.BASE_MODEL_FILE)
    obj_ider = iJR.KAYSER_BIOMASS_IDER
    @info "Orig model" exp size(model) ChU.nzabs_range(model.S)
    
    # Scaling
    scale_factor = 10.0
    @info "Scaling" scale_factor
    model = ChU.well_scaled_model(model, scale_factor)
    
    # fva model
    cfile = joinpath(cache_dir, string("fva_cache_exp_", exp, "sf_", scale_factor, ".bson"))
    if isfile(cfile)
        model = ChU.load_data(cfile)
    else
        model = ChLP.fva_preprocess(model, check_obj = obj_ider, verbose = true)
        ChU.save_data(cfile, model)
    end
    @info "fva model" exp size(model) ChU.nzabs_range(model.S)

    # prepare model
    objidx = ChU.rxnindex(model, obj_ider)
    M, N = size(model)
    exp_growth = Kd.val(:D, exp)
    exp_xi = Kd.val(:xi, exp)
    intake_info = iJR.intake_info(exp)
    ChSS.apply_bound!(model, exp_xi, intake_info; 
        emptyfirst = true)
    fbaout = ChLP.fba(model, objidx)
    fba_growth = ChU.av(model, fbaout, objidx)

    model = ChEP.maxent_ep(model; alpha = Inf, epsconv = 1e-4, 
            maxvar = 1e35, minvar = 1e-35)
end
#     # seed
#     seed_file = iJR.BASE_BETA0_EPOUT_FILE
#     if isfile(seed_file)
#         epout = ChU.load_data(seed_file)
#     else
#         epout = ChEP.maxent_ep(model; alpha = Inf, epsconv = 1e-4, 
#             maxvar = 1e35, minvar = 1e-35)
#         ChU.save_data(seed_file, epout)
#     end

#     # simulation
#     beta_vec = zeros(N)
#     betas = 10.0.^(1:0.1:8)
#     for beta in betas
#         beta_vec[objidx] = beta
#         epout = ChEP.maxent_ep(model; beta_vec, alpha = Inf, damp = 0.9, epsconv = 1e-4, 
#                             maxvar = 1e35, minvar = 1e-35, verbose = true, solution = epout,
#                             maxiter = 1000
#                         )
#         ep_growth = ChU.av(model, epout, objidx)
#         @info "Results" beta exp_growth fba_growth ep_growth
#         println()

#         isnan(ep_growth) && break
#         ep_growth == 0.0 && break
#     end
#     serialize("test.jls", epout)
# end

## -------------------------------------------------------------------
using ProgressMeter
function test_well_scaled_model()
    scale_factors = [1.1, 2.0, 5.0, 10.0, 20.0, 50.0]
    orig_model = ChU.load_data(iJR.BASE_MODEL_FILE)
    obj_ider = iJR.KAYSER_BIOMASS_IDER

    prog = Progress(length(scale_factors))
    for scale_factor in scale_factors

        scl_model = ChU.well_scaled_model(orig_model, scale_factor; 
            verbose = true)

        orig_fbaout = ChLP.fba(orig_model, obj_ider)
        scl_fbaout = ChLP.fba(scl_model, obj_ider)

        for orig_rxn in orig_model.rxns
            scl_rxn = ChU.rxnindex(scl_model, orig_rxn)
            orig_av = ChU.av(orig_model, orig_fbaout, orig_rxn)
            scl_av = ChU.av(scl_model, scl_fbaout, scl_rxn)

            @assert isapprox(orig_av, scl_av; atol = 1e-8)
        end

        next!(prog; showvalues = [
                ("scale_factor        ", scale_factor),
                (": ----------------- ", "ORIGINAL MODEL"),
                ("model size:         ", size(orig_model)),
                ("nzabs_range:        ", ChU.nzabs_range(orig_model.S)),
                ("obj_val:            ", ChU.av(orig_model, orig_fbaout, obj_ider)),
                (": ----------------- ", "SCALED MODEL"),
                ("model size:         ", size(scl_model)),
                ("nzabs_range:        ", ChU.nzabs_range(scl_model.S)),
                ("obj_val:            ", ChU.av(scl_model, scl_fbaout, obj_ider))
            ]
        )

        sleep(3)
    end
    finish!(prog)
end
test_well_scaled_model()