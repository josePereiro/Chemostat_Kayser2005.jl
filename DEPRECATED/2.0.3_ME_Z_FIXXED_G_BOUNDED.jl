let
    method = ME_Z_FIXXED_G_BOUNDED
    objider = iJR.BIOMASS_IDER
    costider = iJR.COST_IDER
    biomass_f = 0.01

    iterator = Kd.val(:D) |> enumerate |> collect 
    @threads for (exp, D) in iterator
        
        thid = threadid()

        # handle cache
        datfile = dat_file(DAT_FILE_PREFFIX; exp, method)
        if isfile(datfile)
            lock(WLOCK) do
                
                @info("Cached loaded (skipping)",
                    exp, D, datfile, threadid()
                ); println()
            end
            continue
        end

        # setup
        model = load_model(exp)
        objidx = ChU.rxnindex(model, objider)
        M, N = size(model)
        exp_growth = Kd.val("D", exp)
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
            

            @info("Finished ", exp, threadid())
            println()
        end
    end
end
