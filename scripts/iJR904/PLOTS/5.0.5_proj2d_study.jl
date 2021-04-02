
## -------------------------------------------------------------------
# proj 2D
let
    method = ME_MAX_POL
    biom_ider = iJR.KAYSER_BIOMASS_IDER

    ps_pool = Dict()
    for exp in EXPS

        for Kd_ider in FLX_IDERS

            # 2D Projection
            p = plot(;title = string("Kayser2005 exp:", exp), 
                xlabel = string(biom_ider), ylabel = string(Kd_ider),
                legend = :left
            )
            proj = DAT[method, :ep, :proj, Kd_ider, exp]
            ChP.plot_projection2D!(p, proj; l = 50)

            # bounds
            lb, ub = DAT[method, :bounds, :flx, Kd_ider, exp]
            hline!(p, [lb]; lw = 3, 
                label = "fva lb",
                color = :blue, ls = :solid
            )
            hline!(p, [ub]; lw = 3,
                label = "fva ub", 
                color = :red, ls = :solid
            )

            # EXPERIMENTAL FLXS
            exp_biom = DAT[method, :Kd, :flx, biom_ider, exp]
            exp_exch = DAT[method, :Kd, :flx, Kd_ider, exp]
            scatter!(p, [exp_biom], [exp_exch]; 
                m = 8, color = :red, label = "exp"
            )
            
            # MAXENT FLXS
            ep_biom = DAT[method, :ep, :flx, biom_ider, exp]
            ep_biom_err = DAT[method, :eperr, :flx, biom_ider, exp]
            ep_exch = DAT[method, :ep, :flx, Kd_ider, exp]
            ep_exch_err = DAT[method, :eperr, :flx, Kd_ider, exp]
            scatter!(p, [ep_biom], [ep_exch]; 
                xerr = [ep_biom_err], yerr = [ep_exch_err],
                m = 8, color = :blue, label = "maxent"
            )

            # mysavefig(p, "polytope"; Kd_ider, exp, method)
            ps_pool[(exp, Kd_ider)] = deepcopy(p)
        end
    end

    # collect 
    for exp in EXPS
        ps = Plots.Plot[ps_pool[(exp, Kd_ider)] for Kd_ider in FLX_IDERS]
        mysavefig(ps, "polytope"; exp, method)
    end

    for Kd_ider in FLX_IDERS
        ps = Plots.Plot[ps_pool[(exp, Kd_ider)] for exp in EXPS]
        mysavefig(ps, "polytope"; Kd_ider, method)
    end
end
