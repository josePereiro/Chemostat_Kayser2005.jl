import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Kayser2005")

import SparseArrays
import Base.Threads: @threads, threadid, SpinLock

# -------------------------------------------------------------------
using Serialization

import Chemostat_Kayser2005
import Chemostat_Kayser2005.Chemostat
import Chemostat_Kayser2005.Chemostat.LP: MathProgBase
const ChK = Chemostat_Kayser2005

const iJR = ChK.iJR904
const Kd = ChK.KayserData # experimental data
const Bd = ChK.BegData    # cost data

const Ch = ChK.Chemostat
const ChU = Ch.Utils
const ChSS = Ch.SteadyState
const ChLP = Ch.LP
const ChEP = Ch.MaxEntEP
const ChSU = Ch.SimulationUtils

using Plots
using Statistics

## -------------------------------------------------------------------
# load data
D = Dict()
for (exp, _) in Kd.val(:D) |> enumerate
    datfile = joinpath(iJR.MODEL_PROCESSED_DATA_DIR, string("ep_dat_exp", exp, ".bson"))
    exp, model, epouts = ChU.load_data(datfile)
    D[exp] = Dict("model" => model, "epouts" => epouts)
end

## -------------------------------------------------------------------
# biomass vs beta
let
    ps = []
    for (exp, _) in Kd.val(:D) |> enumerate

        # p = plot(xlabel = "beta", ylabel = "biomass")
        p = plot(;title = exp, xaxis = nothing, yaxis = nothing)
        model = D[exp]["model"]
        obj_val = iJR.KAYSER_BIOMASS_IDER

        exp_growth = Kd.val(:D, exp)
        fbaout = ChLP.fba(model, obj_val)
        fba_growth = ChU.av(model, fbaout, obj_val)

        for (beta, epout) in D[exp]["epouts"]
            ep_growth = ChU.av(model, epout, obj_val)
            scatter!(p, [beta], [ep_growth]; 
                color = :red, label = "", alpha = 0.5)
        end

        hline!(p, [fba_growth]; color = :blue, 
            ls = :dash, lw = 2, label = "")
        hline!(p, [exp_growth]; color = :black, 
            ls = :dot, lw = 2, label = "")

        push!(ps, p)
    end
    plot(ps...; layout = length(ps), size = (1200, 1200))
end

## -------------------------------------------------------------------
# stoi err vs beta
let
    exp = rand(1:15)
    model = D[exp]["model"]
    M, N = size(model)
    sepouts = sort(collect(D[exp]["epouts"]); by = first)

    # fba
    ex_glc = let
        lmodel = deepcopy(model)
        obj_ider = iJR.KAYSER_BIOMASS_IDER
        ChU.ub!(lmodel, obj_ider, Kd.val(:D, exp))
        fbaout = ChLP.fba(lmodel, obj_ider, iJR.COST_IDER)
        ChU.av(lmodel, fbaout, iJR.GLC_EX_IDER)
    end

    p = plot(title = iJR.PROJ_IDER, 
        xlabel = "beta", ylabel = "log10( |stoi_err| / |ex_glc| )")
    f(x) = log10(x)
    mat = zeros(length(sepouts), M)
    EPS = 1e-8
    for (i, (beta, epout)) in enumerate(sepouts)
        mat[i, :] = (abs.(ChU.stoi_err(model, epout)) ./ abs(ex_glc)) .+ EPS
    end
    plot!(p, first.(sepouts), f.(mat); 
        label = "", color = :black, alpha = 0.1
    )
    plot!(p, first.(sepouts), f.(mean.(eachrow(mat))); 
        label = "mean", color = :red, alpha = 0.8, 
        lw = 5, ls = :dash
    )
    p
end