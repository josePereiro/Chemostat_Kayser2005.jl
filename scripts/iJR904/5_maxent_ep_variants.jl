import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Kayser2005")

@time begin
    import SparseArrays
    import Base.Threads: @threads, threadid

    # -------------------------------------------------------------------
    # Run add https://github.com/josePereiro/Chemostat_Heerden2013.jl in 
    # the Julia Pkg REPL to install the package, then you must activate 
    # the package enviroment (see README)
    import Chemostat_Kayser2005
    const ChK = Chemostat_Kayser2005

    const iJR = ChK.iJR904
    const Kd = ChK.KayserData # experimental data
    const Bd = ChK.BegData    # cost data

    # -------------------------------------------------------------------
    # run add "https://github.com/josePereiro/Chemostat" in the 
    # julia Pkg REPL for installing the package
    import Chemostat
    import Chemostat.LP.MathProgBase
    const Ch = Chemostat
    const ChU = Ch.Utils
    const ChSS = Ch.SteadyState
    const ChLP = Ch.LP
    const ChEP = Ch.MaxEntEP
    const ChSU = Ch.SimulationUtils

    import UtilsJL
    const UJL = UtilsJL

    using Serialization
    # using Statistics
end

## -------------------------------------------------------------------
# globals
const WLOCK = ReentrantLock()
const SIM_GLOBAL_ID = "iJR904_MAXENT_VARIANTS"
const DAT_FILE_PREFFIX =  "maxent_ep_boundle_"

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

## -------------------------------------------------------------------
const HOMO = :HOMO
const BOUNDED = :BOUNDED
const EXPECTED = :EXPECTED

# ## -------------------------------------------------------------------
# let
#     # gradien descent
#     epouts = Dict()
#     x0 = [0.0] 
#     x1 = [1.0]
#     C = [0.5]
#     th = 1e-3
#     epout_seed = nothing
#     target = [10.0]
#     function upfun(x) 
#         @show x
#         sleep(0.1)
#         all(x .< 3) ? 2 .* x : [Inf]
#     end
#     b = ChSU.grad_desc(upfun; x0, x1, th, C, 
#             target, maxiters = 50, verbose = true) |> first
#     @show b
# end

## -------------------------------------------------------------------
# EXPECTED and HOMO
let
    Ds = Kd.val(:D) |> enumerate |> collect
    @threads for (exp, D) in Ds

        # handle cache
        datfile = dat_file(string(DAT_FILE_PREFFIX, EXPECTED); exp)
        dat = isfile(datfile) ? deserialize(datfile) : Dict()

        # setup
        model = base_model(exp)
        objidx = ChU.rxnindex(model, iJR.KAYSER_BIOMASS_IDER)
        M, N = size(model)
        exp_growth = Kd.val(:D, exp)
        growth_ub = ChU.ub(model, iJR.KAYSER_BIOMASS_IDER)
        feasible = exp_growth < growth_ub

        lock(WLOCK) do
            @info("Doing $(EXPECTED)", exp, D, feasible, threadid())
            println()
        end
        !feasible && continue # unfeasible
        
        # gradien descent
        ChU.ub!(model, iJR.KAYSER_BIOMASS_IDER, growth_ub * 1.1) # open a beat helps EP
        epouts = get(dat, :epouts, Dict())
        exp_beta = get(dat, :exp_beta, 0.0) # Must start on zero for HOMO
        x0 = [exp_beta] 
        C0 = 1e5
        x1 = [exp_beta + 100.0]
        th = 1e-2
        target = [exp_growth]
        epout_seed = get(epouts, exp_beta, nothing)
        beta_vec = zeros(size(model, 2))

        upfrec_time = 10 # secunds
        last_uptime = time()
        it = 1

        function upfun(betas)
            beta = first(betas)
            isnan(beta) && return [Inf]

            if haskey(epouts, beta) 
                epout = epouts[beta]
            else
                beta_vec[objidx] = beta
                epout = ChEP.maxent_ep(model; 
                    beta_vec,
                    alpha = Inf,
                    damp = 0.985,
                    epsconv = 1e-4, 
                    verbose = false, 
                    solution = epout_seed,
                    maxiter = 500
                )
            end
            # Error checking
            ep_growth = ChU.av(epout)[objidx]
            isnan(ep_growth) && return [Inf]
            
            # update
            epout_seed = epouts[beta] = epout
            exp_beta = dat[:exp_beta] = beta
            
            update = it == 1 || abs(last_uptime - time()) > upfrec_time || 
                epout.status != ChEP.CONVERGED_STATUS
            if update
                lock(WLOCK) do
                    diff_frac = abs.(exp_growth - ep_growth)/exp_growth
                    @info(
                        "Grad descent... ", 
                        exp, it, 
                        epout.status, epout.iter, 
                        ep_growth, exp_growth, diff_frac, 
                        beta,
                        length(epouts),
                        threadid()
                    ); println()
                    it += 1
                    last_uptime = time()
                end
            end

            # caching
            if rem(it, 10) == 0
                lock(WLOCK) do
                    @info(
                        "Catching... ", 
                        exp, it, 
                        length(epouts),
                        datfile,
                        threadid()
                    ); println()
                end
                serialize(datfile, dat) 
            end

            [ep_growth]
        end

        # maxent
        ret_beta = ChSU.grad_desc(upfun; x0, x1, th, C = [C0], 
            target, maxiters = 150, verbose = false) |> first
            
        lock(WLOCK) do

            isnan(ret_beta) && (@warn("NAN beta", threadid()); println())

            # Storing
            dat[:exp_beta] = exp_beta
            dat[:epouts] = epouts
            dat[:model] = model |> ChU.compressed_model

            # caching
            serialize(datfile, dat)
            INDEX[EXPECTED, :DFILE, exp] = datfile

            ep_growth = ChU.av(epouts[exp_beta])[objidx]
            diff = abs.(exp_growth - ep_growth)
            @info("Finished ",
                it,
                exp, exp_beta, 
                length(epouts),
                ep_growth, exp_growth, diff, 
                threadid()
            ); println()
        end
    end # for (exp, D) in Ds 
end

## -------------------------------------------------------------------
# # BOUNDED
# let
#     biomass_f = 0.01

#     Ds = Hd.val("D") |> enumerate |> collect
#     @threads for (exp, D) in Ds 

#         # handle cache
#         datfile = dat_file(string(DAT_FILE_PREFFIX, BOUNDED); exp)
#         if isfile(datfile)
#             lock(WLOCK) do
#                 INDEX[BOUNDED, :DFILE, exp] = datfile
#                 @info("Cached loaded (skipping)",
#                     exp, D, datfile, threadid()
#                 ); println()
#             end
#             continue
#         end

#         # setup
#         model = base_model(exp)
#         objidx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
#         M, N = size(model)
#         exp_growth = Hd.val("D", exp)
#         ChU.ub!(model, iJR.BIOMASS_IDER, exp_growth * (1.0 + biomass_f))
#         ChU.lb!(model, iJR.BIOMASS_IDER, exp_growth * (1.0 - biomass_f))
#         model = ChLP.fva_preprocess(model, 
#             batchlen = 50,
#             check_obj = iJR.BIOMASS_IDER,
#             verbose = true
#         )

#         lock(WLOCK) do
#             @info("Doing $(BOUNDED)", 
#                 exp, D, threadid()
#             ); println()
#         end

#         # maxent
#         epout = ChEP.maxent_ep(model; alpha = Inf, damp = 0.985, epsconv = 1e-4, 
#                     verbose = false, maxiter = 5000
#                 )
            
#         # storing
#         lock(WLOCK) do
#             # Storing
#             dat = Dict()
#             dat[:exp_beta] = 0.0
#             dat[:epouts] = Dict(0.0 => epout)
#             dat[:model] = model |> ChU.compressed_model

#             # caching
#             serialize(datfile, dat)
#             INDEX[BOUNDED, :DFILE, exp] = datfile

#             @info("Finished ", exp, threadid())
#             println()
#         end
#     end
# end

# ## -------------------------------------------------------------------
# # HOMO
# # It was computed in EXPECTED
# let
#     Ds = Hd.val("D") |> enumerate |> collect
#     @threads for (exp, D) in Ds 

#         lock(WLOCK) do
#             @info("Collecting $(HOMO)", 
#                 exp, D, threadid()
#             ); println()
#         end

#         exp_file = INDEX[EXPECTED, :DFILE, exp]
#         exp_dat = deserialize(exp_file)

#         homo_dat = Dict()
#         homo_dat[:exp_beta] = 0.0
#         epout = exp_dat[:epouts][0.0]  # At beta 0
#         homo_dat[:epouts] = Dict(0.0 => epout)
#         homo_dat[:model] = exp_dat[:model]

#         # save homo
#         homo_file = dat_file(string(DAT_FILE_PREFFIX, HOMO); exp)
#         serialize(homo_file, homo_dat)
#         INDEX[HOMO, :DFILE, exp] = homo_file
            
#     end
# end

# ## -------------------------------------------------------------------
# # FBA
# let
#     # fba
#     # lmodel = deepcopy(model)
#     # ChU.ub!(lmodel, iJR.BIOMASS_IDER, Hd.val("D", exp))
#     # ChLP.fba(lmodel, iJR.BIOMASS_IDER, iJR.COST_IDER)
# end

# ## -------------------------------------------------------------------
# # SAVE INDEX
# ChU.save_data(iJR.MAXENT_VARIANTS_INDEX_FILE, INDEX; verbose = false)