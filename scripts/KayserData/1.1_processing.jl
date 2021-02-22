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
# D vs X
let
    p = plot(; title = "Kayser", 
        xlabel = string("Xv (", Kd.unit(:Xv), ")"), 
        ylabel = string("D (", Kd.unit(:D), ")")
    )
    plot!(p, Kd.val(:D), Kd.val(:Xv); 
        label = "", ls = :dash, lw = 3, color = :black, 
        alpha = 0.6
    )
    scatter!(p, Kd.val(:D), Kd.val(:Xv); 
        label = "", m = 8, color = :black
    )
    mysavefig(p, "X_vs_D") 
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