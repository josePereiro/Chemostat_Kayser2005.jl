## ------------------------------------------------------------------
# ARGS
using ArgParse

set = ArgParseSettings()
@add_arg_table! set begin
    "-w"
        help = "number of workers to use"
        default = "1"
    "--init-clear"
        help = "clear cache before running the simulation"   
        action = :store_true
    "--finish-clear"
        help = "clear cache at the end"   
        action = :store_true
end
parsed_args = parse_args(set)
wcount = parse(Int, parsed_args["w"])
init_clear_flag = parsed_args["init-clear"]
finish_clear_flag = parsed_args["finish-clear"]

## ------------------------------------------------------------------
import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Kayser2005")

using Distributed

## ------------------------------------------------------------------
length(workers()) < wcount && 
    addprocs(wcount; exeflags = "--project")
println("Working in: ", workers())

## ------------------------------------------------------------------
@everywhere begin

    import DrWatson: quickactivate
    quickactivate(@__DIR__, "Chemostat_Kayser2005")

    import Chemostat_Kayser2005: Chemostat
    import Chemostat_Kayser2005.Chemostat.LP: MathProgBase
    const Ch = Chemostat
    const ChU = Chemostat.Utils
    const ChSS = Chemostat.SteadyState
    const ChLP = Chemostat.LP
    const ChSU = Chemostat.SimulationUtils

    import Chemostat_Kayser2005: KayserData, iJO1366
    const iJO = iJO1366
    const Kd = KayserData
    ChU.set_cache_dir(iJO.MODEL_CACHE_DATA_DIR)

end

## ------------------------------------------------------------------
# CLEAR CACHE (WARNING)
if init_clear_flag
    ChU.tagprintln_inmw("CLEARING CACHE ")
    ChU.delete_temp_caches()
    ChU.println_inmw("\n")
end

## ------------------------------------------------------------------
# GLOBAL PARAMS
@everywhere begin

    const sim_params = Dict()
    sim_params[:epochlen] = 50 # This determine how often EP results will be cached

    const epmodel_kwargs = Dict()
    epmodel_kwargs[:alpha] = 1e9

    const epconv_kwargs = Dict()
    epconv_kwargs[:maxiter] = Int(1e4) # The maximum number of iteration before EP to return, even if not converged
    epconv_kwargs[:epsconv] = 1e-7 # The error threshold of convergence
    epconv_kwargs[:maxvar] = 1e10
    epconv_kwargs[:minvar] = 1e-10

    const scaling_params = Dict()
    scaling_params[:scale_base] = 100.0

    params_hash = hash((sim_params, epmodel_kwargs, epconv_kwargs))

end

## ------------------------------------------------------------------
# SIMULATION GLOBAL ID
# This must uniquely identify this simulation version
# It is used to avoid cache collisions
@everywhere sim_global_id = "MAXENT_FBA_EP_v1"

## ------------------------------------------------------------------
# SCALE MODELS
# tagprintln_inmw("SCALING MODELS")
# println_ifmw(" size: ", size(model), " S nzabs_range: ", nzabs_range(model.S), "\n")

## ------------------------------------------------------------------
ChU.tagprintln_inmw("CACHING MODEL")
@everywhere model_hash = (:MODEL, sim_global_id)
model = ChU.load_data(iJO.BASE_MODEL_FILE; verbose = false)
ChU.println_inmw(" size: ", size(model), " S nzabs_range: ", ChU.nzabs_range(model.S), "\n")
ChU.tagprintln_inmw("RESCALING MODEL")
model = ChU.well_scaled_model(model, scaling_params[:scale_base]; verbose = false)
ChU.println_inmw(" size: ", size(model), " S nzabs_range: ", ChU.nzabs_range(model.S), "\n")
ChU.save_cache(model_hash, model; headline = "CACHING MODEL")

## ------------------------------------------------------------------
@everywhere function prepare_model(xi)
    model = ChU.load_cache(model_hash; verbose = false)
    intake_info = iJO.load_base_intake_info()
    # I rescale GLC to reach experimental growth
    intake_info["EX_glc__D_e"]["c"] *= 1.22 
    ChSS.apply_bound!(model, xi, intake_info)
    return model
end

## ------------------------------------------------------------------
# RES IDS
# Collect all the computed results ids for bundling
const chnl = RemoteChannel() do
    Channel{Any}(10)
end
const res_ids = []
const collector = @async while true
    id = take!(chnl)
    push!(res_ids, id)
end

## ------------------------------------------------------------------
# SIMULATION
# Any of the loops can be parallelized by just changing one of the 'map' functions

to_map = Kd.val("D")[1:1] |> enumerate
pmap(to_map) do (Di, D)

    ## SIMULATION PARAMS
    ξs = [Kd.val(:xi, Di)]
    βs = [0.0; ChU.logspace(5.5, 6.5, 50)] 
    
    to_map = Iterators.product(ξs, [βs], Di)
    map(to_map) do (ξ, βs, Di)
        
        ## HASH SEEDS
        # TODO: add PARAMS hash
        sim_hash = (sim_global_id, ξ, Di)

        dat = ChSU.cached_simulation(;
            epochlen = sim_params[:epochlen], 
            # epochlen = 100, # Test
            verbose = true,
            sim_id = sim_hash,
            get_model = function()
                return prepare_model(ξ);
            end,
            objider = iJO.OBJ_IDER, 
            beta_info = [(iJO.OBJ_IDER, βs)],
            costider = iJO.COST_IDER,
            clear_cache = false,
            use_seed = true,
            epmodel_kwargs = epmodel_kwargs,
            epconv_kwargs = epconv_kwargs
        )
        
        ## SAVING DATA
        model = prepare_model(ξ)
        res_id = (:RESULT, sim_hash)
        ChU.save_cache(res_id, (Di, ξ, βs, model, dat); 
            headline = "CATCHING RESULTS\n")
        
        ## PASSING ID TO MASTER
        put!(chnl, res_id)
        
        GC.gc()
        return nothing
    end # map(ξs) do ξ

    return nothing
end; # map(Rd.ststs) do stst

## -------------------------------------------------------------------
# BOUNDLING
sleep(1) # wait for collector to get all ids
bundles = Dict()
for id in res_ids
    Di, ξ, βs, model, dat = ChU.load_cache(id; verbose = false)

    # boundling
    bundle = get!(bundles, Di, ChU.ChstatBundle())

    bundle[ξ, :net] = model
    bundle[ξ, :fba] = dat[:fba]

    for (βi, β) in βs |> enumerate
        bundle[ξ, β, :ep] = dat[(:ep, βi)]
    end
end

## -------------------------------------------------------------------
# SAVING
ChU.save_data(iJO.MAXENT_FBA_EB_BOUNDLES_FILE, bundles)

## ------------------------------------------------------------------
# CLEAR CACHE (WARNING)
if finish_clear_flag
    ChU.tagprintln_inmw("CLEARING CACHE ")
    ChU.delete_temp_caches()
end