import DrWatson
const DW = DrWatson
DW.quickactivate(@__DIR__, "Chemostat_Kayser2005")

# ----------------------------------------------------------------------------
## ARGS
using ArgParse

set = ArgParseSettings()
@add_arg_table! set begin
    "--new-dat"
        help = "ignore disk stored DAT"   
        action = :store_true
    "--skip-me"
        help = "do not recompute MaxEnt part"   
        action = :store_true
    "--skip-lp"
        help = "do not recompute LP part"   
        action = :store_true
end

if isinteractive()
    # Dev values
    new_dat = false
    skip_lp = false
    skip_me = false
else
    parsed_args = parse_args(set)
    new_dat = parsed_args["new-dat"]
    skip_lp = parsed_args["skip-lp"]
    skip_me = parsed_args["skip-me"]
end

# ----------------------------------------------------------------------------
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
# setup container
DAT_FILE = iJR.procdir("dat.bson")
DAT = (!isfile(DAT_FILE) || new_dat) ? ChU.DictTree() : UJL.load_data(DAT_FILE)

# ----------------------------------------------------------------------------------
FLX_IDERS = ["GLC", "CO2", "O2", "AC", "NH4"]
DAT[:FLX_IDERS] = FLX_IDERS

# -------------------------------------------------------------------
# ME methods
const ME_Z_OPEN_G_OPEN        = :ME_Z_OPEN_G_OPEN
const ME_MAX_POL              = :ME_MAX_POL
const ME_MAX_POL_B0           = :ME_MAX_POL_B0
const ME_Z_EXPECTED_G_MOVING  = :ME_Z_EXPECTED_G_MOVING
const ME_Z_EXPECTED_G_BOUNDED = :ME_Z_EXPECTED_G_BOUNDED
const ME_Z_FIXXED_G_BOUNDED   = :ME_Z_FIXXED_G_BOUNDED

# LP methods
const FBA_Z_FIX_MAX_VG_MIN_COST = :FBA_Z_FIX_MAX_VG_MIN_COST
const FBA_Z_FIX_MIN_VG_MIN_COST = :FBA_Z_FIX_MIN_VG_MIN_COST
const FBA_Z_FIX_MAX_VG_MAX_COST = :FBA_Z_FIX_MAX_VG_MAX_COST
const FBA_Z_FIX_MIN_VG_MAX_COST = :FBA_Z_FIX_MIN_VG_MAX_COST
const FBA_Z_VG_FIX_MAX_COST = :FBA_Z_VG_FIX_MAX_COST
const FBA_Z_VG_FIX_MIN_COST = :FBA_Z_VG_FIX_MIN_COST
const FBA_Z_FIX_MAX_COST = :FBA_Z_FIX_MAX_COST
const FBA_Z_FIX_MIN_COST = :FBA_Z_FIX_MIN_COST
const FBA_MAX_Z_MIN_COST = :FBA_MAX_Z_MIN_COST
const FBA_MAX_Z_MAX_COST = :FBA_MAX_Z_MAX_COST


LP_METHODS = [
    FBA_Z_FIX_MAX_VG_MIN_COST, FBA_Z_FIX_MAX_VG_MAX_COST, 
    FBA_Z_FIX_MIN_VG_MIN_COST, FBA_Z_FIX_MIN_VG_MAX_COST,
    FBA_Z_VG_FIX_MAX_COST, FBA_Z_VG_FIX_MIN_COST, FBA_Z_FIX_MAX_COST, 
    FBA_Z_FIX_MIN_COST, FBA_MAX_Z_MIN_COST, FBA_MAX_Z_MAX_COST
]
DAT[:LP_METHODS] = LP_METHODS

ME_METHODS = [
    # ME_Z_OPEN_G_OPEN, 
    ME_MAX_POL,
    ME_MAX_POL_B0,
    # ME_Z_FIXXED_G_BOUNDED, 
    ME_Z_EXPECTED_G_BOUNDED, 
    # ME_Z_EXPECTED_G_MOVING
]
DAT[:ME_METHODS] = ME_METHODS

ALL_METHODS = [LP_METHODS; ME_METHODS]
DAT[:ALL_METHODS] = ALL_METHODS

EXPS = 5:13
DAT[:EXPS] = EXPS

# -------------------------------------------------------------------
ME_INDEX_FILE = iJR.procdir("maxent_ep_index.bson")
ME_INDEX = ChU.load_data(ME_INDEX_FILE; verbose = false);

# -------------------------------------------------------------------
LP_DAT_FILE = iJR.procdir("lp_dat_file.bson")
LP_DAT = ChU.load_data(LP_DAT_FILE; verbose = false);

# ----------------------------------------------------------------------------------
Kd_mets_map = iJR.load_mets_map()
Kd_rxns_map = iJR.load_rxns_map()

# ----------------------------------------------------------------------------------
const DAT_FILE_PREFFIX =  "maxent_ep_dat"

## ----------------------------------------------------------------------------------
# COMMON DAT
let
    max_model = iJR.load_model("max_model"; uncompress = false)
    objider = iJR.BIOMASS_IDER

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
            DAT[:bounds, Kd_ider, exp] = (lb, ub)

        end
    end
end

## ----------------------------------------------------------------------------------
# MAXENT DAT
!skip_me && let 
    WLOCK = ReentrantLock()
    objider = iJR.BIOMASS_IDER

    # util fun
    isexpdep(method) = (method != ME_MAX_POL_B0)
    depks(method, typ, Kd_met, exp) = 
        isexpdep(method) ? (method, typ, Kd_met, exp) : (method, typ, Kd_met)
    function dat_file(;method, exp)
        kwargs = isexpdep(method) ? (;method, exp) : (;method)
        fname = UJL.mysavename(DAT_FILE_PREFFIX, "jls"; kwargs...)
        iJR.procdir(fname)
    end

    # Feed jobs
    nths = nthreads()
    Ch = Channel(nths) do ch
        for method in DAT[:ME_METHODS]
            for exp in DAT[:EXPS]
                put!(ch, (exp, method))
                !isexpdep(method) && break # Do only once
            end
        end
    end

    @threads for thid in 1:nths
        thid = threadid()
        thid == 1 && nths > 1 && continue
        for (exp, method) in Ch

            # ME data
            datfile = dat_file(;method, exp)
            !isfile(datfile) && continue
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
                DAT[depks(method, :flx, "D", exp)...] = ep_biom
                DAT[depks(method, :err, "D", exp)...] = ep_std
            end

            # exchanges
            for Kd_met in FLX_IDERS
                model_exch = Kd_rxns_map[Kd_met]

                # flxs
                ep_av = ChU.av(model, epout, model_exch)
                ep_std = sqrt(ChU.va(model, epout, model_exch))

                # proj 2d
                proj = ChLP.projection2D(model, objider, model_exch; l = 50)
                        
                lock(WLOCK) do
                    DAT[depks(method, :proj, Kd_met, exp)...] = proj
                    DAT[depks(method, :flx, Kd_met, exp)...] = ep_av
                    DAT[depks(method, :err, Kd_met, exp)...] = ep_std
                end
            end

            # additional fluxs
            for (ider, model_iders) in iJR.load_inners_idermap()
                # flxs
                ep_av = ChU.av(model, epout, model_iders[1])
                ep_std = sqrt(ChU.va(model, epout, model_iders[1]))
                if length(model_iders) == 2 # reversible
                    # r = r+ - r-
                    ep_av -= ChU.av(model, epout, model_iders[2])
                    ep_std += sqrt(ChU.va(model, epout, model_iders[2]))
                end

                # proj 2d (fwd only)
                proj = ChLP.projection2D(model, objider, model_iders[1]; l = 50)
                        
                lock(WLOCK) do
                    DAT[depks(method, :proj, ider, exp)...] = proj
                    DAT[depks(method, :flx, ider, exp)...] = ep_av
                    DAT[depks(method, :err, ider, exp)...] = ep_std
                end
            end

        end # for (exp, method)
    end # for thid
end

## ----------------------------------------------------------------------------------
# LP DAT
!skip_lp && let

    objider = iJR.BIOMASS_IDER

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

            # additional fluxs
            for (ider, model_iders) in iJR.load_inners_idermap()
                # flxs
                fba_flx = ChU.av(model, fbaout, model_iders[1])
                if length(model_iders) == 2 # reversible
                    # r = r+ - r-
                    fba_flx -= ChU.av(model, fbaout, model_iders[2])
                end
                        
                DAT[method, :flx, ider, exp] = fba_flx
            end
            
        end

    end # for method
end

## ----------------------------------------------------------------------------------
DAT_FILE = iJR.procdir("dat.bson")
UJL.save_data(DAT_FILE, DAT)