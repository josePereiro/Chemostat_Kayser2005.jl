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
    "--backup-cache"
        help = "create a backup of the cache dir"   
        action = :store_true
    "--ignore-partial-res"
        help = "Ignore the top level caches"   
        action = :store_true
end
parsed_args = parse_args(set)
wcount = parse(Int, parsed_args["w"])
init_clear_flag = parsed_args["init-clear"]
finish_clear_flag = parsed_args["finish-clear"]
backup_cache_flag = parsed_args["backup-cache"]
ignore_partial_res_flag = parsed_args["ignore-partial-res"]

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
end

## ------------------------------------------------------------------
# CACHE DIR
@everywhere begin
    const cache_dir = iJO.MODEL_CACHE_DATA_DIR
    ChU.set_cache_dir(cache_dir)
end

## ------------------------------------------------------------------
# BACKUP CACHES
# TODO: package this
function backup_temp_cache(cache_dir, backup_dir = cache_dir * "_backup"; 
        wtime = 60)
    while true
        files = isdir(cache_dir) ? readdir(cache_dir) : []
        for file in files
            !isdir(backup_dir) && mkpath(backup_dir)
            src_file = joinpath(cache_dir, file)
            dest_file = joinpath(backup_dir, file)
            if !isfile(dest_file) || mtime(src_file) > mtime(dest_file)
                cp(src_file, dest_file; force = true, follow_symlinks = true)
            end
        end
        sleep(wtime)
    end
end
if backup_cache_flag
    ChU.tagprintln_inmw("BACKING UP CACHES \n")
    @async begin
        backup_temp_cache(cache_dir)
    end
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

    # From previous runs
    βintervals = Dict(
        1=>(1.0e6, 4.0e6),
        2=>(2.0e6, 4.5e6),
        3=>(3.5e6, 4.5e6),
        4=>(3.5e6, 4.5e6),
        5=>(3.5e6, 4.5e6),
        6=>(3.5e6, 5.0e6),
        7=>(3e6, 5.06),
        8=>(3e6, 5.5e6),
        9=>(3e6, 5.5e6),
        10=>(3e6, 5.5e6),
        11=>(4e6, 5.5e6),
        12=>(4e6, 5.5e6),
        13=>(4e6, 5.5e6),
        14=>(4e6, 5.5e6),
        15=>(4e6, 5.5e6),
    )

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

to_map = Kd.val("D") |> enumerate
pmap(to_map) do (Di, D)

    ## SIMULATION PARAMS
    ξs = [Kd.val(:xi, Di)]
    β0, β1 = βintervals[Di]
    βs = [0.0; range(β0, β1; length = 25)] 
    # βs = [0.0, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7] 
    
    to_map = Iterators.product(ξs, [βs], Di)
    map(to_map) do (ξ, βs, Di)
        
        # HASH SEEDS
        # TODO: add PARAMS hash
        sim_hash = (sim_global_id, ξ, Di)

        # CHECK CACHE
        res_id = (:RESULT, sim_hash)
        if !ignore_partial_res_flag && isfile(ChU.temp_cache_file(res_id))
            ChU.tagprintln_inmw("RES FILE FOUNDED, CONTINUING ", 
                "\nres_id: ", res_id, "\n")
            ## PASSING ID TO MASTER
            put!(chnl, res_id)
            return nothing
        end

        dat = ChSU.cached_simulation(;
            epochlen = sim_params[:epochlen], 
            # epochlen = 100, # Test
            verbose = true,
            sim_id = sim_hash,
            get_model = function()
                return prepare_model(ξ);
            end,
            objider = iJO.BIOMASS_IDER, 
            beta_info = [(iJO.BIOMASS_IDER, βs)],
            costider = iJO.COST_IDER,
            clear_cache = false,
            use_seed = true,
            epmodel_kwargs = epmodel_kwargs,
            epconv_kwargs = epconv_kwargs
        )
        
        # SAVING DATA
        model = prepare_model(ξ)
        ChU.save_cache(res_id, (Di, ξ, βs, model, dat); 
            headline = "CATCHING RESULTS\n")
        
        # PASSING ID TO MASTER
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

    bundle[ξ, :net] = model |> ChU.compressed_model
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