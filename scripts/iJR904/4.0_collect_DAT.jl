import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Kayser2005")

@time begin
    import SparseArrays
    import Base.Threads: @threads, threadid, SpinLock

    #  ----------------------------------------------------------------------------
    # Run add https://github.com/josePereiro/Chemostat_Kayser2005.jl in the 
    # Julia Pkg REPL to install the package, then you must activate the 
    # package enviroment (see README)
    import Chemostat_Kayser2005
    const ChK = Chemostat_Kayser2005

    const iJR = ChK.iJR904
    const Kd = ChK.KayserData # experimental data
    const Bd = ChK.BegData    # cost data

    #  ----------------------------------------------------------------------------
    # run add "https://github.com/josePereiro/Chemostat" in the 
    # julia Pkg REPL for installing the package
    import Chemostat
    import Chemostat.LP.MathProgBase
    const Ch = Chemostat
    const ChU = Ch.Utils
    const ChLP = Ch.LP

    using Serialization
    
    import UtilsJL
    const UJL = UtilsJL

    using Base.Threads
end

## ----------------------------------------------------------------------------------
DAT = ChU.DictTree();

## ----------------------------------------------------------------------------------
FLX_IDERS = ["GLC", "CO2", "O2", "AC", "NH4"]
DAT[:FLX_IDERS] = FLX_IDERS

## -------------------------------------------------------------------
# ME methods
const ME_Z_OPEN_G_OPEN        = :ME_Z_OPEN_G_OPEN
const ME_MAX_POL              = :ME_MAX_POL
const ME_Z_EXPECTED_G_MOVING  = :ME_Z_EXPECTED_G_MOVING
const ME_Z_EXPECTED_G_BOUNDED = :ME_Z_EXPECTED_G_BOUNDED
const ME_Z_FIXXED_G_BOUNDED   = :ME_Z_FIXXED_G_BOUNDED

# LP methods
const FBA_Z_FIX_MIN_COST      = :FBA_Z_FIX_MIN_COST
const FBA_MAX_BIOM_MIN_COST   = :FBA_MAX_BIOM_MIN_COST
const FBA_Z_FIX_MIN_VG_COST   = :FBA_Z_FIX_MIN_VG_COST
const FBA_Z_VG_FIX_MIN_COST   = :FBA_Z_VG_FIX_MIN_COST

LP_METHODS = [
    FBA_Z_FIX_MIN_COST, FBA_MAX_BIOM_MIN_COST, 
    FBA_Z_FIX_MIN_VG_COST, FBA_Z_VG_FIX_MIN_COST
]
DAT[:LP_METHODS] = LP_METHODS

ME_METHODS = [
    # ME_Z_OPEN_G_OPEN, 
    ME_MAX_POL,
    # ME_Z_FIXXED_G_BOUNDED, 
    # ME_Z_EXPECTED_G_BOUNDED, 
    # ME_Z_EXPECTED_G_MOVING
]
DAT[:ME_METHODS] = ME_METHODS

ALL_METHODS = [LP_METHODS; ME_METHODS]
DAT[:ALL_METHODS] = ALL_METHODS

EXPS = 5:13
DAT[:EXPS] = EXPS

## -------------------------------------------------------------------
ME_INDEX_FILE = iJR.procdir("maxent_ep_index.bson")
ME_INDEX = ChU.load_data(ME_INDEX_FILE; verbose = false);

## -------------------------------------------------------------------
LP_DAT_FILE = iJR.procdir("lp_dat_file.bson")
LP_DAT = ChU.load_data(LP_DAT_FILE; verbose = false);

## ----------------------------------------------------------------------------------
Kd_mets_map = iJR.load_Kd_mets_map()
Kd_rxns_map = iJR.load_Kd_rxns_map()

## ----------------------------------------------------------------------------------
# COMMON DAT
let
    max_model = iJR.load_model("max_model"; uncompress = false)
    objider = iJR.KAYSER_BIOMASS_IDER

    for exp in EXPS
        # exp dat
        Kd_biom = Kd.val("D", exp)
        DAT[:exp, :flx, "D", exp] = Kd_biom
        DAT[:exp, :err, "D", exp] = 0.0

        # bounds
        fva_model = iJR.load_model("fva_models", exp; uncompress = false)
        max_lb, max_ub = ChU.bounds(max_model, objider)
        fva_lb, fva_ub = ChU.bounds(fva_model, objider)
        lb = max(max_lb, fva_lb)
        ub = min(max_ub, fva_ub)
        DAT[:bounds, "D", exp] = (lb, ub)

        for Kd_ider in FLX_IDERS
            model_exch = Kd_rxns_map[Kd_ider]

            # exp dat
            Kd_flx = Kd.uval(Kd_ider, exp)
            DAT[:exp, :flx, Kd_ider, exp] = Kd_flx
            DAT[:exp, :err, Kd_ider, exp] = 0.0

            # bounds
            max_lb, max_ub = ChU.bounds(max_model, model_exch)
            fva_lb, fva_ub = ChU.bounds(fva_model, model_exch)
            lb = max(max_lb, fva_lb)
            ub = min(max_ub, fva_ub)
            DAT[:bounds, model_exch, exp] = (lb, ub)

        end
    end
end

## ----------------------------------------------------------------------------------
# MAXENT DAT
let 
    WLOCK = ReentrantLock()
    objider = iJR.KAYSER_BIOMASS_IDER

    # Feed jobs
    nths = nthreads()
    Ch = Channel(nths) do ch
        for exp in DAT[:EXPS], method in DAT[:ME_METHODS]
            put!(ch, (exp, method))
        end
    end

    @threads for thid in 1:nths
        thid = threadid()
        thid == 1 && nths > 1 && continue
        for (exp, method) in Ch
            
            # ME data
            datfile = ME_INDEX[method, :DFILE, exp]
            datfile == :unfeasible && continue
            dat = deserialize(datfile)
            
            model = dat[:model]
            epouts = dat[:epouts]
            exp_beta = maximum(keys(epouts)) # dat[:exp_beta]
            epout = epouts[exp_beta]
            
            lock(WLOCK) do
                @info("Doing", 
                    exp, method, 
                    length(dat[:epouts]), 
                    epout.iter, thid
                ); println()
            end

            # Biomass
            ep_biom = ChU.av(model, epout, objider)
            ep_std = sqrt(ChU.va(model, epout, objider))
            
            # store
            lock(WLOCK) do
                DAT[method, :flx, "D", exp] = ep_biom
                DAT[method, :err, "D", exp] = ep_std
            end

            for Kd_met in FLX_IDERS
                model_exch = Kd_rxns_map[Kd_met]

                # flxs
                ep_av = ChU.av(model, epout, model_exch)
                ep_std = sqrt(ChU.va(model, epout, model_exch))

                # proj 2d
                proj = ChLP.projection2D(model, objider, model_exch; l = 50)
                        
                lock(WLOCK) do
                    DAT[method, :proj, Kd_met, exp] = proj
                    DAT[method, :flx, Kd_met, exp] = ep_av
                    DAT[method, :err, Kd_met, exp] = ep_std
                end
            end

        end # for (exp, method)
    end # for thid
end

## ----------------------------------------------------------------------------------
# LP DAT
let
    objider = iJR.KAYSER_BIOMASS_IDER

    for method in LP_METHODS
            
        for exp in EXPS

            model = LP_DAT[method, :model, exp]
            fbaout = LP_DAT[method, :fbaout, exp]

            # Biomass
            fba_flx = ChU.av(model, fbaout, objider)
            Kd_flx = Kd.val("D", exp)
            DAT[method, :flx, "D", exp] = fba_flx

            for Kd_ider in FLX_IDERS
                model_ider = Kd_rxns_map[Kd_ider]

                Kd_flx = Kd.uval(Kd_ider, exp)
                fba_flx = ChU.av(model, fbaout, model_ider)
                DAT[method, :flx, Kd_ider, exp] = fba_flx
            end
        end

    end # for method
end

## ----------------------------------------------------------------------------------
DAT_FILE = "dat.bson"
UJL.save_data(DAT_FILE, DAT)