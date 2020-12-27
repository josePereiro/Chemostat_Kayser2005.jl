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
import GR
GR.inline("png")
using Statistics

## -------------------------------------------------------------------
# load data
D = Dict()
for (exp, _) in Kd.val(:D) |> enumerate
    datfile = joinpath(iJR.MODEL_PROCESSED_DATA_DIR, string("ep_dat_exp", exp, ".bson"))
    exp, model, epouts, fbaout = ChU.load_data(datfile)
    D[exp] = Dict("model" => model, "epouts" => epouts, "fbaout" => fbaout)
end

# -------------------------------------------------------------------
# Meta
sname = "2.1" # source name
fig_path(name) = joinpath(iJR.MODEL_FIGURES_DIR, name)

## -------------------------------------------------------------------
# biomass vs beta
let
    ps = []
    for (exp, _) in Kd.val(:D) |> enumerate

        p = plot(;xlabel = "beta", ylabel = "biomass", 
            title = string(iJR.PROJ_IDER, " exp_", exp))
        model = D[exp]["model"]
        obj_val = iJR.KAYSER_BIOMASS_IDER

        exp_growth = Kd.val(:D, exp)
        fbaout = D[exp]["fbaout"]
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

        # saving
        fname = fig_path(string(sname, "_exp_", exp, "_biomass_vs_beta.png"))
        savefig(p, fname)
    end
end

## -------------------------------------------------------------------
# stoi err vs beta
function _plot_stoi_err_vs_beta(exp)
    
    model = D[exp]["model"]
    M, N = size(model)
    sepouts = sort(collect(D[exp]["epouts"]); by = first)

    fbaout = D[exp]["fbaout"]
    ex_glc = ChU.av(model, fbaout, iJR.GLC_EX_IDER)

    p = plot(title = string(iJR.PROJ_IDER, " exp_", exp), 
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
    # saving
    fname = fig_path(string(sname, "_exp_", exp, "_stoi_err_vs_beta.png"))
    savefig(p, fname)
    p
end
for (exp, _) in Kd.val(:D) |> enumerate
    _plot_stoi_err_vs_beta(rand(1:15))
end

## -------------------------------------------------------------------
function get_closest_beta(exp, obj_ider)
    
    exp_val = Kd.val(:D, exp)
    model = D[exp]["model"]
    sepouts = sort(collect(D[exp]["epouts"]); by = first)

    closest_beta = nothing
    closest_diff = Inf
    for (beta, epout) in sepouts
        ep_val = ChU.av(model, epout, obj_ider)
        diff = abs(ep_val - exp_val)
        if closest_diff > diff
            closest_diff = diff
            closest_beta = beta
        end
    end
    closest_beta
end

## -------------------------------------------------------------------
# biomass correlation
let
    obj_ider = iJR.KAYSER_BIOMASS_IDER

    fba_plot = plot(;xlabel = "exp biomass", ylabel = "model biomass", 
        title = string(iJR.PROJ_IDER, " FBA"))
    ep_plot = plot(;xlabel = "exp biomass", ylabel = "model biomass", 
        title = string(iJR.PROJ_IDER, " EP"))

    min_, max_ = Inf, -Inf
    for (exp, _) in Kd.val(:D) |> enumerate
        
        @show exp

        model = D[exp]["model"]
        M, N = size(model)
        sepouts = sort(collect(D[exp]["epouts"]); by = first)
        exp_beta = get_closest_beta(exp, obj_ider)
        epout = D[exp]["epouts"][exp_beta]
        fbaout = D[exp]["fbaout"]
        
        # norm
        fba_val = ChU.av(model, fbaout, obj_ider)
        ep_val = ChU.av(model, epout, obj_ider)
        ep_err = sqrt(ChU.va(model, epout, obj_ider))
        exp_val = Kd.val(:D, exp , nothing)
        isnothing(exp_val) && continue
        @show exp_val fba_val ep_val
        scatter!(fba_plot, [exp_val], [fba_val]; 
            label = "", color = :black, alpha = 0.5, ms = 5)

        ep_err = sqrt(ChU.va(model, epout, obj_ider))
        scatter!(ep_plot, [exp_val], [ep_val]; yerr = [ep_err],
            label = "", color = :black, alpha = 0.5, ms = 5)

        
        min_ = minimum([min_, exp_val, fba_val, ep_val])
        max_ = maximum([max_, exp_val, fba_val, ep_val])
    end

    for plt in [fba_plot, ep_plot]
        plot!(plt, [min_,max_], [min_,max_]; label = "", lw = 2, alpha = 0.4, ls = :dash, color = :black)
    end
    p = plot([fba_plot, ep_plot]...; layout = 2)
    fname = fig_path(string(sname, "_biomass_correlation.png"))
    savefig(p, fname)
    nothing
end

## -------------------------------------------------------------------
# full exchs correlation
let
    # TODO: package this
    Kd_msd_mets = ["GLC", "CO2", "O2", "AC", "NH4"]
    met_map = iJR.load_mets_map()
    exch_map = iJR.load_exch_met_map()
    obj_ider = iJR.KAYSER_BIOMASS_IDER

    fba_exch_plot = plot(;xlabel = "exp exch / |ex_GLC|", ylabel = "model exch / |ex_GLC|", 
        title = string(iJR.PROJ_IDER, " FBA"))
    ep_exch_plot = plot(;xlabel = "exp exch / |ex_GLC|", ylabel = "model exch / |ex_GLC|", 
        title = string(iJR.PROJ_IDER, " EP"))

    min_, max_ = Inf, -Inf
    for (exp, _) in Kd.val(:D) |> enumerate
        
        @show exp

        model = D[exp]["model"]
        M, N = size(model)
        sepouts = sort(collect(D[exp]["epouts"]); by = first)
        exp_beta = get_closest_beta(exp, obj_ider)
        epout = D[exp]["epouts"][exp_beta]
        fbaout = D[exp]["fbaout"]
        
        # norm
        ex_glc = Kd.uval("GLC", exp, nothing)
        isnothing(ex_glc) && continue
        @show ex_glc

        for Kd_met in Kd_msd_mets
            iJR_met = met_map[Kd_met]
            iJR_ex = exch_map[iJR_met]
            fba_val = ChU.av(model, fbaout, iJR_ex)  / abs(ex_glc)
            ep_val = ChU.av(model, epout, iJR_ex)  / abs(ex_glc)
            ep_err = sqrt(ChU.va(model, epout, iJR_ex)) / abs(ex_glc)
            exp_val = -Kd.uval(Kd_met, exp , nothing) / abs(ex_glc)
            isnothing(exp_val) && continue
            
            @show exp_val fba_val ep_val
            scatter!(fba_exch_plot, [exp_val], [fba_val]; 
                label = "", color = :black, alpha = 0.5)

            scatter!(ep_exch_plot, [exp_val], [ep_val]; 
                yerr = [ep_err],
                label = "", color = :black, alpha = 0.5)
            
            min_ = minimum([min_, exp_val, fba_val, ep_val])
            max_ = maximum([max_, exp_val, fba_val, ep_val])
        end
    end

    @show min_ max_
    for plt in [fba_exch_plot, ep_exch_plot]
        plot!(plt, [min_,max_], [min_,max_]; label = "", lw = 2, alpha = 0.4, ls = :dash, color = :black)
    end
    p = plot([fba_exch_plot, ep_exch_plot]...; layout = 2)
    fname = fig_path(string(sname, "_exchanges_correlation.png"))
    savefig(p, fname)
    nothing
end

## -------------------------------------------------------------------
# full conc correlation
let
    # TODO: package this
    Kd_msd_mets = ["GLC", "AC", "NH4"]
    met_map = iJR.load_mets_map()
    exch_map = iJR.load_exch_met_map()
    obj_ider = iJR.KAYSER_BIOMASS_IDER

    fba_conc_plot = plot(;xlabel = "exp conc / cGLC", ylabel = "model conc / cGLC", 
        title = string(iJR.PROJ_IDER, " FBA"))
    ep_conc_plot = plot(;xlabel = "exp conc / cGLC", ylabel = "model conc / cGLC", 
        title = string(iJR.PROJ_IDER, " EP"))

    min_, max_ = Inf, -Inf
    for (exp, _) in Kd.val(:D) |> enumerate
        
        @show exp

        model = D[exp]["model"]
        M, N = size(model)
        sepouts = sort(collect(D[exp]["epouts"]); by = first)
        exp_beta = get_closest_beta(exp, obj_ider)
        epout = D[exp]["epouts"][exp_beta]
        fbaout = D[exp]["fbaout"]
        exp_xi = Kd.val("xi", exp, nothing)
        isnothing(exp_xi) && continue
        # norm
        cglc = Kd.cval("GLC", exp, nothing)
        isnothing(cglc) && continue
        @show cglc

        for Kd_met in Kd_msd_mets
            iJR_met = met_map[Kd_met]
            iJR_ex = exch_map[iJR_met]
            fba_flx = ChU.av(model, fbaout, iJR_ex)
            ep_flx = ChU.av(model, epout, iJR_ex)
            ep_flxerr = sqrt(ChU.va(model, epout, iJR_ex)) 
            
            
            @show fba_flx ep_flx

            # s = c + u*xi
            exp_c = Kd.cval(Kd_met, exp , 0.0)
            exp_s = Kd.sval(Kd_met, exp , nothing)
            isnothing(exp_s) && continue
            exp_s /= cglc
            fba_s = max(0.0, (exp_c + fba_flx * exp_xi)) / cglc
            ep_s = max(0.0, exp_c + ep_flx * exp_xi)/ cglc
            ep_serr = max(0.0, ep_flxerr * exp_xi) / cglc 

            scatter!(fba_conc_plot, [exp_s], [fba_s]; 
                label = "", color = :black, alpha = 0.5)

            scatter!(ep_conc_plot, [exp_s], [ep_s]; 
                yerr = [ep_serr],
                label = "", color = :black, alpha = 0.5)
            
            min_ = minimum([min_, exp_s, fba_s, ep_s])
            max_ = maximum([max_, exp_s, fba_s, ep_s])
        end
    end

    @show min_ max_
    for plt in [fba_conc_plot, ep_conc_plot]
        plot!(plt, [min_,max_], [min_,max_]; label = "", lw = 2, alpha = 0.4, ls = :dash, color = :black)
    end
    p = plot([fba_conc_plot, ep_conc_plot]...; layout = 2)
    fname = fig_path(string(sname, "_conc_correlation.png"))
    savefig(p, fname)
    nothing
end