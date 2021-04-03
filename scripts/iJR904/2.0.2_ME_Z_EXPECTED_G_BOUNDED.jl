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
            model = load_model(exp)
            objidx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
            M, N = size(model)
            exp_growth = Kd.val(:D, exp)
            growth_ub = ChU.ub(model, iJR.BIOMASS_IDER)
            feasible = exp_growth < growth_ub
            biom_lb, biom_ub = ChU.bounds(model, iJR.BIOMASS_IDER)
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
            ChU.ub!(model, iJR.BIOMASS_IDER, growth_ub * 1.1) # open a beat helps EP

            lock(WLOCK) do
                nzabs_range = ChU.nzabs_range(model.S)
                @info("Starting... ", 
                    exp, method,
                    size(model), nzabs_range, 
                    feasible,
                    threadid()
                ); println()
            end
            !feasible && continue

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
                dat[:exp_beta] = maximum(keys(epouts))
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
    objider = iJR.BIOMASS_IDER

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
