## -------------------------------------------------------------------
# fva bounds
let
   
    ps = Plots.Plot[]
    for ider in FLX_IDERS
        p = plot(title = ider, xlabel = "replica", ylabel = "flx")
        xticks =  (EXPS, string.(EXPS))
        
        Kd_vals = DAT[:Kd, :flx, ider, EXPS]
        plot!(p, EXPS, Kd_vals; 
            label = "exp", color = :black, alpha = 0.8, lw = 3, xticks)

        for method in ALL_METHODS
            color = method_colors[method]    

            bounds = DAT[method, :bounds, :flx, ider, EXPS]
            lb, ub = first.(bounds), last.(bounds)
            plot!(p, EXPS, [lb ub]; 
                fillrange = ub,
                fillalpha = 0.2,
                label = "", color = :gray
            )
            
            ep_vals = DAT[method, :ep, :flx, ider, EXPS]
            ep_errs = DAT[method, :eperr, :flx, ider, EXPS]
            plot!(p, EXPS, ep_vals; 
                label = string(method), color, alpha = 0.5, 
                fillrange = [ep_vals .- ep_errs ep_vals .+ ep_errs],
                lw = 1, xticks
            )
        end
        mysavefig(p, "bound_study"; ider)
        push!(ps, p)
    end
    mysavefig(ps, "bound_study")
    
end