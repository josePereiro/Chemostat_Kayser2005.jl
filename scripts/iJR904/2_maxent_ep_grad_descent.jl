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
    for (exp, D) in Kd.val(:D) |> enumerate

        model = ChU.load_data(iJR.BASE_MODEL_FILE)
        obj_ider = iJR.KAYSER_BIOMASS_IDER
        @info "Orig model" exp size(model) ChU.nzabs_range(model.S)

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

        # simulation
        datfile = joinpath(iJR.MODEL_PROCESSED_DATA_DIR, string("ep_dat_exp", exp, ".bson"))
        epouts = isfile(datfile) ? ChU.load_data(datfile).epouts : Dict()
        beta_vec = zeros(N)

        # log approach
        last_beta = maximum(keys(epouts))
        epout_seed = isempty(epouts) ? nothing : epouts[last_beta]
        betas = [0.0; 10.0.^(3:0.05:8)]
        nan_beta = last_beta

        for approach in [:log_approach, :linear_approach]
            
            @info "Starting" exp approach

            for beta in betas

                nan_beta = beta
                haskey(epouts, beta) && continue

                beta_vec[objidx] = beta
                epout = ChEP.maxent_ep(model; 
                                    beta_vec, alpha = Inf, damp = 0.9, epsconv = 1e-4, 
                                    maxvar = 1e50, minvar = 1e-50, verbose = true, solution = epout_seed,
                                    maxiter = 1000
                                )
                # info
                ep_growth = ChU.av(model, epout, objidx)
                @info "Results" exp D beta exp_growth fba_growth ep_growth
                println()

                # out conditions
                isnan(ep_growth) && break
                ep_growth == 0.0 && break

                # updating
                epout_seed = epout
                epouts[beta] = epout
            end

            # lineal approach approach
            last_beta = maximum(keys(epouts))
            betas = range(last_beta, 2 * nan_beta; length = 1000)
        end

        # saving
        ChU.save_data(datfile, (;exp, model, epouts))
    
    end # for (exp, D)
end

