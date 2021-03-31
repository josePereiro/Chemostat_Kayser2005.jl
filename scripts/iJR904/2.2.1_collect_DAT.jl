DAT = ChU.DictTree()
let 

    WLOCK = ReentrantLock()

    # CACHE
    DATfile = joinpath(iJR.MODEL_PROCESSED_DATA_DIR, "2.1_DAT.jls")
    if isfile(DATfile) 
        global DAT = deserialize(DATfile) 
        @info("DAT CACHE LOADED")
        return
    end
    DAT[:EXPS] = []

    objider = iJR.KAYSER_BIOMASS_IDER
    DAT[:CONC_IDERS] = CONC_IDERS
    DAT[:FLX_IDERS] = FLX_IDERS

    # Find exps
    for exp in Kd.EXPS
        ok = false
        for method in ALL_METHODS
            ok = haskey(INDEX, method, :DFILE, exp) &&
                INDEX[method, :DFILE, exp] != :unfeasible
            !ok && break
        end
        !ok && continue
        push!(DAT[:EXPS], exp)
    end
    max_model = iJR.load_model("max_model"; uncompress = false)

    # Feed jobs
    nths = nthreads()
    Ch = Channel(nths) do ch
        for exp in DAT[:EXPS], method in ALL_METHODS
            put!(ch, (exp, method))
        end
    end


    @threads for thid in 1:nths
        thid = threadid()
        thid == 1 && nths > 1 && continue
        for (exp, method) in Ch
                
            !haskey(INDEX, method, :DFILE, exp) && continue
            datfile = INDEX[method, :DFILE, exp]
            datfile == :unfeasible && continue
            dat = deserialize(datfile)
            
            model = dat[:model]
            objidx = ChU.rxnindex(model, objider)
            epouts = dat[:epouts]
            exp_beta = maximum(keys(epouts)) # dat[:exp_beta]
            epout = epouts[exp_beta]
            exp_xi = Kd.val(:xi, exp)
            fva_model = iJR.load_model("fva_models", exp; uncompress = false)
            
            lock(WLOCK) do
                @info("Doing", 
                    exp, method, 
                    length(dat[:epouts]), 
                    epout.iter, thid
                ); println()
            end

            # Biomass
            ep_biom = ChU.av(model, epout, objidx)
            ep_std = sqrt(ChU.va(model, epout, objidx))
            Kd_biom = Kd.val("D", exp)
            max_lb, max_ub = ChU.bounds(max_model, objidx)
            fva_lb, fva_ub = ChU.bounds(fva_model, objidx)
            lb = max(max_lb, fva_lb)
            ub = min(max_ub, fva_ub)
            
            # store
            lock(WLOCK) do
                DAT[method, :ep   , :flx, objider, exp] = ep_biom
                DAT[method, :eperr, :flx, objider, exp] = ep_std
                DAT[method, :Kd   , :flx, objider, exp] = Kd_biom
                DAT[:Kd   , :flx, objider, exp] = Kd_biom
                DAT[method, :bounds, :flx, objider, exp] = (lb, ub)
            end

            # fluxes
            for Kd_met in FLX_IDERS

                    model_met = Kd_mets_map[Kd_met]
                    model_exch = Kd_rxns_map[Kd_met]
                    model_exchi = ChU.rxnindex(model, model_exch)

                    ep_av = ChU.av(model, epout, model_exchi)
                    ep_std = sqrt(ChU.va(model, epout, model_exchi))
                    Kd_flx = Kd.val("u$Kd_met", exp)
                    proj = ChLP.projection2D(model, objider, model_exchi; l = 50)
                    
                    max_lb, max_ub = ChU.bounds(max_model, Kd_rxns_map[Kd_met])
                    fva_lb, fva_ub = ChU.bounds(fva_model, Kd_rxns_map[Kd_met])
                    lb = max(max_lb, fva_lb)
                    ub = min(max_ub, fva_ub)
                            
                    lock(WLOCK) do
                        DAT[method, :Kd, :flx, Kd_met, exp] = Kd_flx
                        DAT[:Kd, :flx, Kd_met, exp] = Kd_flx
                        DAT[method, :ep, :proj, Kd_met, exp] = proj
                        DAT[method, :ep, :flx, Kd_met, exp] = ep_av
                        DAT[method, :eperr, :flx, Kd_met, exp] = ep_std
                        DAT[method, :bounds, :flx, Kd_met, exp] = (lb, ub)
                    end

            end

            # mets
            for Kd_met in CONC_IDERS

                ep_std = DAT[method, :eperr, :flx, Kd_met, exp] 
                ep_av = DAT[method, :ep, :flx, Kd_met, exp]
                # conc (s = c + u*xi)
                c = Kd.val("c$Kd_met", exp, 0.0)
                ep_conc = max(c + ep_av * exp_xi, 0.0)
                Kd_conc = Kd.val("s$Kd_met", exp)

                lock(WLOCK) do
                    DAT[method, :Kd, :conc, Kd_met, exp] = Kd_conc
                    DAT[:Kd, :conc, Kd_met, exp] = Kd_conc
                    DAT[method, :ep, :conc, Kd_met, exp] = ep_conc
                    DAT[method, :eperr, :conc, Kd_met, exp] = ep_std * exp_xi
                end
            end

        end # for (exp, method)
    end # for thid

    # saving
    DAT[:EXPS] |> unique! |> sort!
    serialize(DATfile, DAT)
end

# -------------------------------------------------------------------
# Inter project comunication
let
    CORR_DAT = isfile(iJR.CORR_DAT_FILE) ? ChU.load_data(iJR.CORR_DAT_FILE) : Dict()
    CORR_DAT[:MAXENT_EP] = DAT
    ChU.save_data(iJR.CORR_DAT_FILE, CORR_DAT)
end
