import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Kayser2005")

@time begin
    import SparseArrays
    import Base.Threads: @threads, threadid, SpinLock

    using Serialization

    import Chemostat_Kayser2005
    const ChK = Chemostat_Kayser2005

    const iJR = ChK.iJR904
    const Kd = ChK.KayserData # experimental data
    const Bd = ChK.BegData    # cost data

    import Chemostat
    import Chemostat.LP: MathProgBase
    const Ch = ChK.Chemostat
    const ChU = Ch.Utils
    const ChSS = Ch.SteadyState
    const ChLP = Ch.LP
    const ChEP = Ch.MaxEntEP
    const ChSU = Ch.SimulationUtils

    import UtilsJL
    const UJL = UtilsJL
    using Serialization
    using Base.Threads
    UJL.set_cache_dir(iJR.MODEL_CACHE_DIR)
end

## -------------------------------------------------------------------
# globals
const WLOCK = ReentrantLock()
const SIM_GLOBAL_ID = "iJR904_MAXENT_VARIANTS"
const DAT_FILE_PREFFIX =  "maxent_ep_dat_"

const INDEX = UJL.DictTree()
function dat_file(name; kwargs...)
    fname = UJL.mysavename(name, "jls"; kwargs...)
    joinpath(iJR.MODEL_PROCESSED_DATA_DIR, fname)
end

## -------------------------------------------------------------------
function base_model(exp)
    BASE_MODELS = ChU.load_data(iJR.BASE_MODELS_FILE; verbose = false);
    model_dict = BASE_MODELS["fva_models"][exp]
    ChU.MetNet(;model_dict...) |> ChU.uncompressed_model
end

# -------------------------------------------------------------------
# METHOD VARIANTS
const ME_Z_OPEN_G_OPEN          = :ME_Z_OPEN_G_OPEN           # Do not use extra constraints
const ME_Z_EXPECTED_G_BOUNDED   = :ME_Z_EXPECTED_G_BOUNDED    # Match ME and Dy biom average and constraint av_ug
const ME_Z_EXPECTED_G_MOVING    = :ME_Z_EXPECTED_G_MOVING     # 
const ME_Z_FIXXED_G_BOUNDED     = :ME_Z_FIXXED_G_BOUNDED      # Fix biom around observed

## -------------------------------------------------------------------
# ME_Z_EXPECTED_G_MOVING
let
    reset_cache = false
    method = ME_Z_EXPECTED_G_MOVING

    # Feed jobs
    Ch = Channel(1) do ch
        cGLCs = Kd.val("cGLC")
        for (exp, cGLC)  in enumerate(cGLCs)
            put!(ch, (exp, cGLC))
        end
    end

    # thid = 1 # Test
    # exp, cGLC = 3, Kd.val("cGLC", 3) # Test
    @threads for thid in 1:nthreads()
        for (exp, cGLC) in Ch

            ## -------------------------------------------------------------------
            # handle cache
            datfile = dat_file(string(DAT_FILE_PREFFIX, method); exp)
            if isfile(datfile)
                lock(WLOCK) do
                    INDEX[method, :DFILE, exp] = datfile
                    @info("Cached loaded (skipping)",
                        exp, cGLC,datfile, thid
                    )
                    println()
                end
                continue
            end
            
            ## -------------------------------------------------------------------
            # SetUp
            model_cache_id = (:MODEL0_CACHE, exp)
            reset_cache && UJL.delete_cache(model_cache_id) # uncomment to reset
            model =  UJL.load_cache(model_cache_id; verbose = true) do
                BASE_MODELS = ChU.load_data(iJR.BASE_MODELS_FILE; verbose = false);
                model_dict = BASE_MODELS["fva_models"][exp]
                model0 = ChU.MetNet(;model_dict...) |> ChU.uncompressed_model

                scale_factor = 0.5
                for rxni in eachindex(model0.rxns)
                    lb, ub = ChU.bounds(model0, rxni)
                    Δ = ub - lb
                    ChU.lb!(model0, rxni, lb - Δ * scale_factor)
                    ChU.ub!(model0, rxni, ub + Δ * scale_factor)
                end

                return ChLP.fva_preprocess(model0; 
                    verbose = false, check_obj = iJR.KAYSER_BIOMASS_IDER
                )   
            end
            objidx = ChU.rxnindex(model, iJR.KAYSER_BIOMASS_IDER)
            exp_growth = Kd.val("D", exp)
            biom_lb, biom_ub = ChU.bounds(model, iJR.KAYSER_BIOMASS_IDER)
            if biom_ub < exp_growth
                lock(WLOCK) do
                    INDEX[method, :DFILE, exp] = :unfeasible
                    @info("Not feasible (skipping)", 
                        biom_ub ,exp_growth, 
                        thid
                    ); println()
                end
                continue
            end
            cgD_X = -Kd.cval(:GLC, exp) * Kd.val(:D, exp) / Kd.val(:Xv, exp)
            exglc_L = ChU.lb(model, iJR.GLC_EX_IDER)
            exglc_qta = abs(exglc_L * 0.005)
            expβ = 0.0
            vg_avPME = 0.0
            epouts_cid = (:EPOUTS_CACHE, exp)
            reset_cache && UJL.delete_cache(epouts_cid) 
            epouts = ChU.load_cache(epouts_cid, Dict(); verbose = false)
            epout_cid = (:EPOUT_CACHE, exp)
            reset_cache && UJL.delete_cache(epout_cid) 
            epout = ChU.load_cache(epout_cid; verbose = false)

            for movround in 1:500

                ## -------------------------------------------------------------------
                # GRAD DESCEND
                x0 = expβ
                x1 = 10.0
                maxΔ = max(expβ * 0.05, 1e3)
                gd_th = 1e-2
                target = exp_growth
                beta_vec = zeros(size(model, 2))
        
                upfrec_time = 50 # secunds
                last_uptime = time()
                gd_it = 1
        
                ## -------------------------------------------------------------------
                function upfun(beta)
        
                    beta_vec[objidx] = beta
                    epouts[beta] = epout = ChEP.maxent_ep(model; 
                        beta_vec,
                        alpha = Inf,
                        maxiter = 2000, 
                        epsconv = 1e-3, 
                        verbose = false, 
                        solution = epout
                    )

                    # ep_growth = ChU.av(epout)[objidx]
                    ep_growth = ChU.av(model, epout, iJR.KAYSER_BIOMASS_IDER)
                    vg_avPME = ChU.av(model, epout, iJR.GLC_EX_IDER)
        
                    update = gd_it == 1 || abs(last_uptime - time()) > upfrec_time || 
                        epout.status != ChEP.CONVERGED_STATUS
        
                    update && lock(WLOCK) do
                        diff = abs.(exp_growth - ep_growth)
                        @info(
                            "Grad descent... ", 
                            exp, gd_it, 
                            epout.status, epout.iter, 
                            ep_growth, exp_growth, diff, 
                            (biom_lb, biom_ub),
                            beta, thid
                        ); println()
                        last_uptime = time()
        
                        try;
                            ChU.save_cache(epouts_cid, epouts; verbose = false)
                            ChU.save_cache(epout_cid, epout; verbose = false)
                        catch err
                            @warn("ERROR SAVING IGNORED", err)
                        end
                    end
                    
                    gd_it += 1
                    return ep_growth
                end
        
                ## -------------------------------------------------------------------
                # FIND BETA
                expβ = UJL.grad_desc(upfun; x0, x1, th = gd_th, maxΔ, 
                    target, 
                    maxiters = 5000, 
                    verbose = false
                )
        
                ## -------------------------------------------------------------------
                # MOVE V_UB
                Δstep = 0.5
                exglc_lb, exglc_ub = ChU.bounds(model, iJR.GLC_EX_IDER)

                # vg_avPME = exglc_lb * 0.8
                # lb is the uptake limit
                dist = cgD_X - vg_avPME
                Δexglc_lb = sign(dist) * max(exglc_qta, abs(dist * Δstep))
                exglc_lb = min(exglc_ub,
                    min(cgD_X, 
                        max(exglc_L, exglc_lb + Δexglc_lb)
                    )
                )
                ChU.lb!(model, iJR.GLC_EX_IDER, exglc_lb)
                
                ## -------------------------------------------------------------------
                # INFO AND CONV
                conv = cgD_X <= vg_avPME && epout.status == ChEP.CONVERGED_STATUS
                
                lock(WLOCK) do
                    @info("Round Done", 
                        movround, conv,
                        dist, exglc_qta, Δexglc_lb,
                        (vg_avPME, cgD_X), 
                        exglc_ub, exglc_lb,  exglc_L, 
                        thid
                    ); println()
                end
                conv && break
        
            end #  for movround in 1:1000

            ## -------------------------------------------------------------------
            lock(WLOCK) do

                # Storing
                dat = Dict()
                dat[:exp_beta] = expβ
                dat[:epouts] = epouts
                dat[:model] = model |> ChU.compressed_model

                # caching
                serialize(datfile, dat)
                INDEX[method, :DFILE, exp] = datfile

                ep_growth = ChU.av(epouts[expβ])[objidx]
                diff = abs.(exp_growth - ep_growth)
                @info("Finished ",
                    exp, expβ, 
                    length(epouts),
                    ep_growth, exp_growth, diff, 
                    thid
                ); println()
            end

        end # for (exp, cGLC) in Ch
    end # for thid in 1:nthreads()
end

## -------------------------------------------------------------------
# ME_Z_EXPECTED_G_BOUNDED
let
    # global setup
    method = ME_Z_EXPECTED_G_BOUNDED

    # Feed jobs
    Ch = Channel(1) do ch
        iterator = Kd.val(:D) |> enumerate |> collect 
        for (exp, D) in iterator
            put!(ch, (exp, D))
        end
    end

    @threads for thid in 1:nthreads()
        for (exp, D) in Ch

            # prepare model
            model = base_model(exp)
            objidx = ChU.rxnindex(model, iJR.KAYSER_BIOMASS_IDER)
            M, N = size(model)
            exp_growth = Kd.val(:D, exp)
            growth_ub = ChU.ub(model, iJR.KAYSER_BIOMASS_IDER)
            feasible = exp_growth < growth_ub
            biom_lb, biom_ub = ChU.bounds(model, iJR.KAYSER_BIOMASS_IDER)
            if biom_ub < exp_growth
                lock(WLOCK) do
                    INDEX[method, :DFILE, exp] = :unfeasible
                    @info("Not feasible (skipping)", 
                        exp, method,
                        biom_ub ,exp_growth, 
                        thid
                    ); println()
                end
                continue
            end
            ChU.ub!(model, iJR.KAYSER_BIOMASS_IDER, growth_ub * 1.1) # open a beat helps EP

            lock(WLOCK) do
                nzabs_range = ChU.nzabs_range(model.S)
                @info("Starting... ", 
                    exp, method,
                    size(model), nzabs_range, 
                    feasible,
                    threadid()
                ); println()
            end
            # !feasible && continue

            # simulation
            datfile = dat_file(DAT_FILE_PREFFIX; exp, method)
            dat = isfile(datfile) ? deserialize(datfile) : Dict()
            epouts = get!(dat, :epouts, Dict())
            init_len = length(epouts)
            beta_vec = zeros(N)
            approach_status = get!(dat, :approach_status, :running)
            if approach_status == :finished 
                lock(WLOCK) do
                    INDEX[method, :DFILE, exp] = datfile
                    @show("FINISHED", 
                        exp, method,
                        approach_status, 
                        datfile, 
                        thid
                    ); println()
                end
                continue
            end
            convth = 0.05
            
            # log approach
            epout_seed = isempty(epouts) ? nothing : epouts[maximum(keys(epouts))]
            betas = [0.0; 10.0.^(3:0.05:15)]
            nan_beta = first(betas)

            for approach in [:log_approach, :linear_approach]
                
                lock(WLOCK) do
                    @info("Starting", 
                        exp, approach, 
                        length(epouts),
                        threadid()
                    ); println()
                end

                for beta in betas

                    nan_beta = beta
                    haskey(epouts, beta) && continue

                    beta_vec[objidx] = beta
                    epout = nothing
                    try
                        epout = ChEP.maxent_ep(model; 
                            beta_vec, alpha = Inf, damp = 0.9, epsconv = 1e-4, 
                            maxvar = 1e50, minvar = 1e-50, verbose = false, solution = epout_seed,
                            maxiter = 1000
                        )
                    catch err; end

                    # info
                    ep_growth = isnothing(epout) ? 0.0 : ChU.av(model, epout, objidx)
                    lock(WLOCK) do
                        @info("Results", exp, method, 
                            beta, 
                            exp_growth, growth_ub, ep_growth, 
                            length(epouts),
                            threadid()
                        ); println()
                    end

                    # error conditions
                    fail = isnothing(epout) || isnan(ep_growth) || ep_growth == 0.0 
                    fail && break

                    # updating
                    epout_seed = epouts[beta] = epout

                    # convergence
                    converr = abs(ep_growth - exp_growth)/exp_growth
                    conv = converr < convth || ep_growth > exp_growth 
                    conv && break

                end # for betas

                # Catching
                # update = init_len != length(epouts)
                update = true
                update && lock(WLOCK) do
                    serialize(datfile, dat)
                    @info("Catching", 
                        exp, method,  
                        length(epouts),
                        basename(datfile),
                        threadid()
                    ); println()
                end

                # lineal approach
                last_beta = maximum(keys(epouts))
                step = abs(last_beta - nan_beta) / 100.0
                iszero(step) && break
                betas = range(last_beta, 1.0e15; step)

            end # for approach

            # saving
            lock(WLOCK) do
                INDEX[method, :DFILE, exp] = datfile
                dat[:approach_status] = :finished
                dat[:model] = model |> ChU.compressed_model
                serialize(datfile, dat)
                @info("Finished", 
                    exp, method,  
                    length(epouts),
                    threadid()
                ); println()
            end
       
        end # for (exp, D)
    end
end


## -------------------------------------------------------------------
# Further convergence
let
    method = ME_Z_EXPECTED_G_BOUNDED
    objider = iJR.KAYSER_BIOMASS_IDER

    iterator = Kd.val(:D) |> enumerate |> collect 
    @threads for (exp, D) in iterator

        datfile = INDEX[method, :DFILE, exp]
        datfile == :unfeasible && continue
        dat = deserialize(datfile)
        model, epouts = ChU.uncompressed_model(dat[:model]) , dat[:epouts]
        
        exp_growth = Kd.val(:D, exp)
        exp_beta = maximum(keys(epouts))
        exp_epout = epouts[exp_beta]

        lock(WLOCK) do
            @info("Converging...", 
                exp, method,
                exp_beta, exp_epout.status, 
                threadid()
            ); println()
        end
        converg_status = get!(dat, :converg_status, :undone)
        converg_status == :done && continue

        new_epout = nothing
        if exp_epout.status == :unconverged
            try;
                objidx = ChU.rxnindex(model, objider)
                beta_vec = zeros(size(model, 2)); 
                beta_vec[objidx] = exp_beta
                new_epout = ChEP.maxent_ep(model; 
                    beta_vec, alpha = Inf, damp = 0.98, epsconv = 1e-4, 
                    maxvar = 1e50, minvar = 1e-50, verbose = false, 
                    solution = exp_epout, maxiter = 5000
                )
            catch err; @warn("ERROR", err); println() end

            ep_growth = isnothing(new_epout) ? 0.0 : ChU.av(model, new_epout, objider)
            fail = isnan(ep_growth) || ep_growth == 0.0 
            epouts[exp_beta] = fail ? exp_epout : new_epout
        end
        
        # Saving
        lock(WLOCK) do
            @info("Saving...", 
                exp, method, 
                exp_beta, 
                exp_epout.status,
                threadid()
            ); println()
        end
        dat[:converg_status] = :done
        serialize(datfile, dat)
    end
end

## -------------------------------------------------------------------
# ME_Z_FIXXED_G_BOUNDED
let
    method = ME_Z_FIXXED_G_BOUNDED
    objider = iJR.KAYSER_BIOMASS_IDER
    costider = iJR.COST_IDER
    biomass_f = 0.01

    iterator = Kd.val(:D) |> enumerate |> collect 
    @threads for (exp, D) in iterator
        
        thid = threadid()

        # handle cache
        datfile = dat_file(DAT_FILE_PREFFIX; exp, method)
        if isfile(datfile)
            lock(WLOCK) do
                INDEX[method, :DFILE, exp] = datfile
                @info("Cached loaded (skipping)",
                    exp, D, datfile, threadid()
                ); println()
            end
            continue
        end

        # setup
        model = base_model(exp)
        objidx = ChU.rxnindex(model, objider)
        M, N = size(model)
        exp_growth = Kd.val("D", exp)
        biom_lb, biom_ub = ChU.bounds(model, iJR.KAYSER_BIOMASS_IDER)
        if biom_ub < exp_growth
            lock(WLOCK) do
                INDEX[method, :DFILE, exp] = :unfeasible
                @info("Not feasible (skipping)", 
                    exp, method,
                    biom_ub ,exp_growth, 
                    thid
                ); println()
            end
            continue
        end
        
        fbaout = ChLP.fba(model, objider, costider)
        fba_growth = ChU.av(model, fbaout, objider)
        ub_growth = min(fba_growth, exp_growth)
        ChU.ub!(model, objider, ub_growth * (1.0 + biomass_f))
        ChU.lb!(model, objider, ub_growth * (1.0 - biomass_f))
        model = ChLP.fva_preprocess(model, 
            check_obj = objider,
            verbose = false
        )

        lock(WLOCK) do
            @info("Doing... ", 
                exp, method, 
                D, threadid()
            ); println()
        end

        # maxent
        epout = ChEP.maxent_ep(model; 
            alpha = Inf, damp = 0.985, epsconv = 1e-4, 
            verbose = false, maxiter = 5000
        )
            
        # storing
        lock(WLOCK) do
            # Storing
            dat = Dict()
            dat[:epouts] = Dict(0.0 => epout)
            dat[:model] = model |> ChU.compressed_model

            # caching
            serialize(datfile, dat)
            INDEX[method, :DFILE, exp] = datfile

            @info("Finished ", exp, threadid())
            println()
        end
    end
end

## -------------------------------------------------------------------
# ME_Z_OPEN_G_OPEN
let
    method = ME_Z_OPEN_G_OPEN
    objider = iJR.KAYSER_BIOMASS_IDER

    iterator = Kd.val(:D) |> enumerate |> collect 
    @threads for (exp, D) in iterator

        thid = threadid()

        # handle cache
        datfile = dat_file(DAT_FILE_PREFFIX; exp, method)
        if isfile(datfile)
            lock(WLOCK) do
                INDEX[method, :DFILE, exp] = datfile
                @info("Cached loaded (skipping)",
                    exp, D, datfile, threadid()
                ); println()
            end
            continue
        end

        # setup
        model = base_model(exp)
        objidx = ChU.rxnindex(model, objider)

        exp_growth = Kd.val("D", exp)
        biom_lb, biom_ub = ChU.bounds(model, iJR.KAYSER_BIOMASS_IDER)
        if biom_ub < exp_growth
            lock(WLOCK) do
                INDEX[method, :DFILE, exp] = :unfeasible
                @info("Not feasible (skipping)", 
                    exp, method,
                    biom_ub ,exp_growth, 
                    thid
                ); println()
            end
            continue
        end

        lock(WLOCK) do
            @info("Doing... ", 
                exp, method, 
                D, threadid()
            ); println()
        end

        # maxent
        epout = ChEP.maxent_ep(model; 
            alpha = Inf, damp = 0.985, epsconv = 1e-4, 
            verbose = false, maxiter = 5000
        )
            
        # storing
        lock(WLOCK) do
            # Storing
            dat = Dict()
            dat[:epouts] = Dict(0.0 => epout)
            dat[:model] = model |> ChU.compressed_model

            # caching
            serialize(datfile, dat)
            INDEX[method, :DFILE, exp] = datfile

            @info("Finished ", exp, threadid())
            println()
        end
    end
end

## -------------------------------------------------------------------
# save index
ChU.save_data(iJR.MAXENT_VARIANTS_INDEX_FILE, INDEX; verbose = false)