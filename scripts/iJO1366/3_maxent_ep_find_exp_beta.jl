import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Kayser2005")

## ------------------------------------------------------------------

import SparseArrays
import Chemostat_Kayser2005: Chemostat
import Chemostat_Kayser2005.Chemostat.LP: MathProgBase
const Ch = Chemostat
const ChU = Chemostat.Utils
const ChSS = Chemostat.SteadyState
const ChLP = Chemostat.LP
const ChSU = Chemostat.SimulationUtils
const ChEP = Chemostat.MaxEntEP
import Chemostat_Kayser2005: KayserData, iJO1366
const iJO = iJO1366
const Kd = KayserData

## ------------------------------------------------------------------
# LOAD DATA
bundles = ChU.load_data(iJO.MAXENT_FBA_EB_BOUNDLES_FILE)

## ------------------------------------------------------------------
let
    Di = 1
    xi = Kd.val("xi", Di)
    bundle = bundles[Di]
    model = bundle[xi, :net]
    model = ChU.uncompressed_model(model)
    M, N = size(model)
    iders = [iJO.BIOMASS_IDER]
    idxs = [ChU.rxnindex(model, ider) for ider in iders]
    fbaout = ChLP.fba(model, iJO.BIOMASS_IDER, iJO.COST_IDER)
    fba_cost_val = ChU.av(model, fbaout, iJO.COST_IDER)
    ChU.bounds!(model, iJO.COST_IDER, fba_cost_val * 0.9, fba_cost_val * 1.1)
    ChU.summary(model, iJO.COST_IDER)
    target = [ChU.av(model, fbaout, iJO.BIOMASS_IDER)]
    x0 = [0.0] # Initial beta
    x1 = [1e1] # first step
    C = [1e2]
    th = 1e-3
    epouts = Dict()
    epouts[x0] = bundle[xi, 0.0, :ep]
    epout_seed = epouts[x0]
    @info string("Seed objval: ", ChU.av(model, epout_seed, iJO.BIOMASS_IDER))
    @info string("target: ", target)
    beta_vec = zeros(N)
    running = true
    t = @async while running
        file = joinpath(@__DIR__, "epouts.bson")
        ChU.save_data(file, epouts; verbose = false)
        sleep(30)
    end
    expb = ChSU.grad_desc(;x0, x1, th, C, target, maxiters = 5000) do betas
        if haskey(epouts, betas) 
            epout = epouts[betas]
        else
            beta_vec[idxs] = betas
            epsconv = 1e-4
            alpha = 1e8
            verbose = false
            solution = epout_seed
            epout = ChEP.maxent_ep(model; alpha,  beta_vec, epsconv, verbose, solution)
            epouts[betas] = epout
        end
        epout_seed = epout
        r = ChU.av(epout)[idxs]
    end
    running = false
    wait(t)
end
## ------------------------------------------------------------------


