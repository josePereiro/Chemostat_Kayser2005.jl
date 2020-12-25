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
cproj = Base.active_project()
length(workers()) < wcount && 
    addprocs(wcount; exeflags = "--project=$cproj")
println("Working in: ", workers())

## ------------------------------------------------------------------
@everywhere begin

    # import DrWatson: quickactivate
    # quickactivate(@__DIR__, "Chemostat_Kayser2005")

    import Chemostat_Kayser2005
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
    epmodel_kwargs[:alpha] = Inf

    const epconv_kwargs = Dict()
    epconv_kwargs[:maxiter] = Int(1e5) # The maximum number of iteration before EP to return, even if not converged
    epconv_kwargs[:epsconv] = 1e-4 # The error threshold of convergence
    epconv_kwargs[:maxvar] = 1e45
    epconv_kwargs[:minvar] = 1e-45

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
ChU.tagprintln_inmw("CACHING MODEL")
@everywhere model_hash = (:MODEL, sim_global_id)
model = ChU.load_data(iJO.BASE_MODEL_FILE; verbose = false)
model = ChU.clampfields!(model, [:b, :lb, :ub]; abs_max = iJO.ABS_MAX_BOUND, zeroth = iJO.ZEROTH)
ChU.println_inmw(" size: ", size(model), " S nzabs_range: ", ChU.nzabs_range(model.S), "\n")
# ChU.tagprintln_inmw("RESCALING MODEL")
# model = ChU.well_scaled_model(model, scaling_params[:scale_base]; verbose = false)
# ChU.println_inmw(" size: ", size(model), " S nzabs_range: ", ChU.nzabs_range(model.S), "\n")
ChU.save_cache(model_hash, model; headline = "CACHING MODEL")

## ------------------------------------------------------------------
@everywhere function prepare_model(xi)
    model = ChU.load_cache(model_hash; verbose = false)
    intake_info = iJO.load_base_intake_info()
    # I rescale GLC to reach experimental growth
    intake_info["EX_glc__D_e"]["c"] *= 1.25 
    ChSS.apply_bound!(model, xi, intake_info)
    # for cost to be around fba minimal value
    fbaout = ChLP.fba(model, iJO.BIOMASS_IDER, iJO.COST_IDER)
    fba_cost_val = ChU.av(model, fbaout, iJO.COST_IDER)
    ChU.bounds!(model, iJO.COST_IDER, 0.0, fba_cost_val * 1.1)
    return model
end

## ------------------------------------------------------------------
let # Test fba
    for (Di, exp_objval) in enumerate(Kd.val("D"))
        xi = Kd.val("xi", Di)
        model = prepare_model(xi)
        fbaout = ChLP.fba(model, iJO.BIOMASS_IDER, iJO.COST_IDER)
        fba_objval = ChU.av(model, fbaout, iJO.BIOMASS_IDER)
        @info "FBA Test" Di xi fba_objval exp_objval
        println()
    end
end

## ------------------------------------------------------------------
# RES IDS
# Collect all the computed results ids for bundling
const chnl = RemoteChannel() do
    Channel{Any}(length(workers()))
end
const res_ids = []
const collector = @async begin
    while true
        !isopen(chnl) && isempty(chnl) && break
        push!(res_ids, take!(chnl))
    end
end

## ------------------------------------------------------------------
# SIMULATION
# Any of the loops can be parallelized by just changing one of the 'map' functions
model = nothing; GC.gc()
to_map = Kd.val("D") |> enumerate
pmap(to_map) do (Di, D)

    ## SIMULATION PARAMS
    ξs = [Kd.val(:xi, Di)]
    βs = [0.0] 
    
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
            verbose = true,
            sim_id = sim_hash,
            get_model = () -> prepare_model(ξ),
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
close(chnl)
wait(collector) # wait for collector to get all ids

## -------------------------------------------------------------------
# BOUNDLING
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