let
    method = ME_Z_OPEN_G_OPEN
    objider = iJR.BIOMASS_IDER

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
        model = load_model(exp)
        objidx = ChU.rxnindex(model, objider)

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