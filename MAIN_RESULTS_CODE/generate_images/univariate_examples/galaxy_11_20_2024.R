# setwd("/Users/hb23255/Dropbox/DSS_MIX/DSS_MIX_06_2024/")
setwd("/Users/hbolfarine/Dropbox/DSS_MIX/DSS_MIX_06_2024/")
source("source/dcpossum_dens_comp.R")
source("source/dcpossum_sim_data_mix.R")
source("source/dcpossum_plots.R")
source("source/dcpossum_unc.R")
source("source/func_pred_laplace_temp.R")
source("source/dcpossum_clust.R")
# source("source/dcpossum_clust_compar.R")


#----------------------------------------------##----------------------------------------------#
#----------------------------------------------##----------------------------------------------#

# Approved - needs to change the distances of the points related to the classification
# plot two pictures, K = 3, K = 4

#----------------------------------------------##----------------------------------------------#
#----------------------------------------------##----------------------------------------------#

y.data.app = data_sim_func("galaxy")
set.seed(1800) 
# k0 = 0.1 - galaxy
# DPM.galaxy = dcpossum.DPM(y.data.app, kmax = 10, quant.sample = 1000)
DPM.galaxy = dcpossum.DPM.dir(y.data.app, kmax = 10, quant.sample = 1000, k0 = 1/10, pred.f = TRUE)
K_star = 4

# DPM_comp_galaxy = plot.possum.uni(DPM.galaxy[[1]], kmax = 10, y.data = y.data.app, sel.K = FALSE, y.lim = c(-1.2,0.5))
DPM_comp_galaxy = plot.possum.uni(DPM.galaxy[[1]], kmax = 10, sel.K = FALSE, y.lim = c(-1.1,0.4))
DPM_comp_galaxy

possum.DPM.galaxy = possum.unc.quant.values(DPM.galaxy, K.sel = K_star)

plots.possum.exemp.DPM = plot.possum.quant(possum.DPM.galaxy, K.sel = K_star, scale.plot = FALSE, model = "DPM",
                                       y.data.app, index.possum = FALSE, index.pred = "galaxy")
plots.possum.exemp.DPM$dens.summ

# DPM_comp_galaxy + plots.possum.exemp$dens.summ

possum_clust = dc.possum.clust.uni(DPM.galaxy, y.data.app, K.sel = K_star, km = TRUE)

dat.clust.DPM = process_clustering_uni(DPM.galaxy[[6]], possum_clust, K_star, y.data.app)

plot.galaxy.DPM = create_custom_plot(dat.clust.DPM, plots.possum.exemp.DPM$dens.summ,
                                    k_star = K_star, y.c = -0.012, y.min = -0.1, y.max = 1, y.k = -0.03, text_plot = "Velocities of Galaxies")



# comp_DPM = plot_posterior_components(table(DPM.galaxy[[3]]))

# galaxy_DPM = DPM_comp_galaxy + plot.galaxy.DPM 
# 
# # plot for k = 4
# ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/galaxy_DPM.pdf", 
#        plot = galaxy_DPM, width = 12, height = 3.7)

# DPM.galaxy[[2]]
# temp = densityMclust(DPM.galaxy[[6]],G = 4,modelNames ="E")
# temp$parameters
# temp$classification
# 
# plot(DPM.galaxy[[6]], temp$classification)

#----------------------------------------------##----------------------------------------------#
#----------------------------------------------##----------------------------------------------#

y.data.app = data_sim_func("galaxy")
set.seed(1800)

MFM.galaxy = dcpossum.MFM(y.data.app, kmax = 10, quant.sample = 1000, pred.f = TRUE)

K_star = 4
MFM_comp = plot.possum.uni(MFM.galaxy[[1]], kmax = 10, sel.K = FALSE,y.lim = c(-1.3,0.5))
MFM_comp

possum.MFM.galaxy = possum.unc.quant.values(MFM.galaxy, K.sel = K_star)

plots.possum.exemp.MFM = plot.possum.quant(possum.MFM.galaxy, K.sel = K_star, y.data.app, scale.plot = FALSE,
                                           model = "MFM",index.possum = FALSE, index.pred = "galaxy")

# MFM_comp + plots.possum.exemp.MFM$dens.summ

# # # Plot - number of components from the methods

possum_clust = dc.possum.clust.uni(MFM.galaxy, y.data.app, K.sel = K_star, km = TRUE)

dat.clust.MFM = process_clustering_uni(MFM.galaxy[[6]], possum_clust, K_star, y.data.app)

# y = -0.07, -0.15
plot.galaxy.MFM = create_custom_plot(dat.clust.MFM, plots.possum.exemp.MFM$dens.summ, k_star = K_star, y.c = -0.012, y.k = -0.03,
                                     text_plot = "Velocities of Galaxies")

# MFM_post_clust = plot_posterior_components(table(MFM.galaxy[[3]]))

# galaxy_MFM = MFM_comp + plot.galaxy.MFM

# plot for k = 4
# ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/galaxy_MFM.pdf", 
#        plot = galaxy_MFM, width = 12, height = 3.7)

#----------------------------------------------##----------------------------------------------#
#----------------------------------------------##----------------------------------------------#

set.seed(1822) 
y.data.app = data_sim_func("galaxy")
BP.galaxy = dcpossum.BP(y.data.app, kmax = 10, BP.run = F, quant.sample = 1000, data.bp = "galaxy", pred.f = TRUE, scale.plot = c(10,40))

K_star = 3

# temp = densityMclust(BP.galaxy[[6]], G = 3)
# print(temp)

BP_comp_galaxy = plot.possum.uni(BP.galaxy[[1]], kmax = 10, sel.K = FALSE, y.lim = c(-1.2,0.5))
BP_comp_galaxy

possum.BP.galaxy = possum.unc.quant.values(BP.galaxy, K.sel = K_star)

# setwd("/Users/hb23255/")
plots.possum.exemp = plot.possum.quant(possum.BP.galaxy, K.sel = K_star, scale.plot = FALSE, 
                                       y.data.app, index.possum = FALSE, index.pred = "galaxy", model = "BP")
plots.possum.exemp$dens.summ

# BP_comp_galaxy + plots.possum.exemp$dens.summ

possum_clust = dc.possum.clust.uni(BP.galaxy, y.data.app, K.sel = K_star, km = TRUE)

dat.clust.BP = process_clustering_uni(BP.galaxy[[6]], possum_clust, K_star, y.data.app)

# dat.clust.BP.3  = dat.clust.BP
# dens.3 = plots.possum.exemp$dens.summ
# y = -0.07, -0.15
plot.galaxy.BP.3 = create_custom_plot(dat.clust.BP, plots.possum.exemp$dens.summ, 
                                    k_star = K_star, y.c = -0.012, y.k = -0.025, text_plot = "Velocities of Galaxies")

plot.galaxy.BP.4 = create_custom_plot(dat.clust.BP, plots.possum.exemp$dens.summ, 
                                      k_star = K_star, y.c = -0.012, y.k = -0.025, text_plot = "Velocities of Galaxies")

plot.galaxy.BP.final.3 = ggdraw() +
  draw_plot(plot.galaxy.BP.3) +  # Your main plot
  draw_plot(plots.possum.exemp$legend, x = 0.1, y = 0.35, width = 0.35, height = 0.95)

plot.galaxy.BP.final.4 = ggdraw() +
  draw_plot(plot.galaxy.BP.4) +  # Your main plot
  draw_plot(plots.possum.exemp$legend, x = 0.1, y = 0.35, width = 0.35, height = 0.95)

# plot for k = 3
# ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/paper_test_overleaf/galaxy_BP_ksel3.pdf", plot = plot.galaxy.BP.final.3, width = 6, height = 4.5)
# 
# # plot for k = 4
# ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/paper_test_overleaf/galaxy_BP_ksel4.pdf", plot = plot.galaxy.BP.final.4, width = 6, height = 4.5)

# plot_posterior_components(table(DPM.galaxy[[3]]))

# Both plots
plot.galaxy.BP.both = plot.galaxy.BP.final.3 + plot.galaxy.BP.final.4 + plot_annotation(tag_levels = list(c("(a)", "(b)"))) 
ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/sn-article-template/plot.galaxy.BP.both.pdf", plot = plot.galaxy.BP.both, width = 12, height = 4.5)

comp.all = plot.postcomp.uni(DPM.galaxy,MFM.galaxy,BP.galaxy,K.true = FALSE)

# BP_model_eval = BP_comp_galaxy + comp.all + plot_annotation(tag_levels = list(c("(a)", "(b)")))

# plot for k = 4
ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/galaxy_possum_post_clust.pdf", 
       plot = BP_model_eval, width = 12, height = 3.5)

#----------------------------------------------##----------------------------------------------#
#----------------------------------------------##----------------------------------------------#

# Plot 2

comp_galaxy = plot.possum.comp.uni(BP.galaxy,DPM.galaxy,MFM.galaxy, kmax = 10)
comp_galaxy + comp.all

p_comp_galaxy = (comp_galaxy + comp.all + plot_layout(widths = c(2, 1))) + plot_annotation(tag_levels = list(c("a)", "b)")))
# ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/sim_plot_discr.pdf", width = 11, height = 3.5)
# 

# plot.galaxy.DPM + plot.galaxy.MFM

p_comp_pred_galaxy = ((comp_galaxy + comp.all) + plot_layout(widths = c(2, 1))) / (plot.galaxy.DPM + plot.galaxy.MFM) + 
  plot_annotation(tag_levels = list(c("(a)", "(b)","(c)","(d)"))) +
  plot_layout(heights = c(1.8, 2.5))

p.legend <- ggdraw() +
  draw_plot(p_comp_pred_galaxy) +  # Your main plot
  draw_plot(plots.possum.exemp.DPM$legend, x = 0.05, y = 0.31, width = 0.25, height = 0.25) +
  draw_plot(plots.possum.exemp.MFM$legend, x = 0.53, y = 0.31, width = 0.25, height = 0.25)
   
ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/galaxy_mult_all.pdf", 
       plot = p.legend, width = 12.5, height = 7.5)
 


