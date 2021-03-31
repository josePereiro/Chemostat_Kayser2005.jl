import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Kayser2005")

@time begin
    
    import Chemostat_Kayser2005
    const ChK = Chemostat_Kayser2005
    const Kd = ChK.KayserData    

    using UtilsJL
    const UJL = UtilsJL

    using Plots
end

## -------------------------------------------------------------------
fileid = "1.1"
function mysavefig(p, pname; params...) 
    fname = UJL.mysavefig(p, string(fileid, "_", pname), Kd.KAYSER_FIGURES_DIR; params...)
    @info "Plotting" fname
end

## -------------------------------------------------------------------
# BALANCE
let
    ps = Plots.Plot[]
    for met in [:GLC, :NH4]
        p = plot(; title = string("Balance: ", met), 
            xlabel = "feed", 
            ylabel = "exch + drain" 
        )

        exps = 1:13 

        feed = Kd.cval.(met, exps) .* Kd.val.(:D, exps)
        exch = Kd.uval.(met, exps) .* Kd.val.(:Xv, exps) .|> abs
        drain = Kd.sval.(met, exps) .* Kd.val.(:D, exps) .|> abs

        scatter!(p, feed, exch .+ drain; 
            label = "", m = 8, color = :black
        )
        vals = [feed; exch .+ drain] |> sort
        plot!(p, vals, vals;
            label = "", ls = :dash, lw = 3, alpha = 0.6, color = :black
        )
        push!(ps, p)
    end
    mysavefig(ps, "Balances") 
end

## -------------------------------------------------------------------
function plot_ids(id1, id2)
    xs, ys = Kd.val(id1), Kd.val(id2)
    ml = min(length(xs), length(ys))
    xs, ys = xs[1:ml], ys[1:ml]
    sids = sortperm(xs)
    ys = ys[sids]
    p = scatter(;xlabel = string(id1), ylabel = string(id2))
    scatter!(p, xs, ys; label = "", color = :black, m = 8)
    plot!(p, xs, ys; label = "", color = :black, lw = 3, alpha = 0.5)
end

## -------------------------------------------------------------------
# corrs
let
    ids = [:D, :sAC, :sGLC, :sNH4, :Xv, :uAC, :uGLC, :uNH4]
    for i in 1:length(ids)
        for j in (i + 1):length(ids)
            id1 = ids[i]
            id2 = ids[j]
            p = plot_ids(id1, id2)
            mysavefig(p, "ids_plot"; id1, id2) 
        end
    end
end 

## -------------------------------------------------------------------
# gasses extrapolation
# the gases exchanges are non reported 
# for the experiments [14, 15]...
# I computed from an extrapolation with the CTR and OTR.
let
    ps = Plots.Plot[]
    for (uid, tid) in [("uCO2", "CTR"), ("uO2", "OTR")]
        u = Kd.TABLE2[uid][3:end] # mmol/ g hr
        tr = Kd.TABLE1[tid][3:end] # g/ L hr
        idxs = 1:min(length(u), length(tr))
        xs = tr[idxs]
        ys = u[idxs]
        p = scatter(xs, ys; 
            label = "reported", 
            xlabel = string(tid, " (g/ L hr)"), 
            ylabel = string(uid, " (mmol/ g hr)"), 
            m = 8, alpha = 0.8, legend = :left
        )
        m = sum(ys ./ xs)/ length(xs)
        plot!(p, xs, m .* xs; 
            label = "", lw = 3, ls = :dash, 
            alpha = 0.6
        )
        
        # extrapolation
        xs = tr
        ys = u
        eys = m .* xs
        scatter!(p, xs, eys; 
            label = "extrap", alpha = 0.5, m = 8
        )
        push!(ps, p)
    end
    mysavefig(ps, "tr_vs_u")
end

## -------------------------------------------------------------------
# exchage for exp 14 and 15 wans't reported
# I compute it using 0 = uX + (c - s)D
let
    ps = Plots.Plot[]
    for id in [:AC, :GLC, :NH4]
        com_us = []
        rep_us = []
        p = plot(;xlabel = "exp", ylabel = string(id))
        for exp in Kd.EXPS
            c = Kd.cval(id, exp, 0.0)
            s = Kd.sval(id, exp)
            D = Kd.val(:D, exp)
            Xv = Kd.val(:Xv, exp)
            cu = -(c - s) * D / Xv
            eu = Kd.uval(id, exp)
            push!(com_us, cu)
            push!(rep_us, eu)
        end
        scatter!(p, exps, com_us; 
            label = "computed", m = 8, alpha = 0.5
        )
        scatter!(p, exps, rep_us; 
            label = "reported", m = 8, alpha = 0.5
        )
        push!(ps, p)
    end

    mysavefig(ps, "exchanges_computed")

end