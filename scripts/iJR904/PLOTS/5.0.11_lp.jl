# # -----------------------------------------------------------------------------------------------
# # Flx correlations
# let
#     ps = Plots.Plot[]
#     for method in ALL_METHODS
#         p = plot(;title = string(method), xlabel = "exp flx", ylabel = "model flx")
#         for Kd_ider in FLX_IDERS
#             model_ider = Kd_rxns_map[Kd_ider]
#             for exp in EXPS
#                 color = ider_colors[Kd_ider]
#                 # every body is possitive here
#                 Kd_flx = abs(Kd.uval(Kd_ider, exp))

#                 model = LPDAT[method, :model, exp]
#                 fbaout = LPDAT[method, :fbaout, exp]
                
#                 fba_flx = abs(ChU.av(model, fbaout, model_ider))
#                 DAT[method, :Kd, :flx, Kd_ider, exp] = Kd_flx
#                 DAT[method, :lp, :flx, Kd_ider, exp] = fba_flx

#                 scatter!(p, [Kd_flx], [fba_flx]; 
#                     ms = 8, color, label = ""
#                 )
#             end
#         end
#         xs = DAT[method, [:Kd, :lp], :flx, FLX_IDERS, EXPS] |> sort
#         plot!(p, xs, xs; label = "", ls = :dash, 
#             alpha = 0.9, color = :black, lw = 3
#         )
#         push!(ps, p)
#     end
    
#     pname = "flx_tot_corr"
#     mysavefig(ps, pname)

# end

## -------------------------------------------------------------------
# # # leyends
# # # TODO fix this...
# # let
# #     for (title, colors) in [
# #             ("exp", exp_colors), 
# #             ("iders", ider_colors),
# #             ("method", method_colors)
# #         ]
# #     p = plot(; framestyle = :none)
# #         scolors = sort(collect(colors); by = (p) -> string(first(p)))
# #         for (id, color) in scolors
# #             scatter!(p, [0], [0];
# #                 thickness_scaling = 1,
# #                 color, ms = 8, label = string(id),
# #                 legendfontsize=10, 
# #                 # size = [300, 900],
# #                 # legend = :left
# #             )
# #         end
# #         mysavefig(p, "$(title)_color_legend")
# #     end
# # end