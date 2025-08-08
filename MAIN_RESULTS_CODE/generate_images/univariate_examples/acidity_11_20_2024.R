setwd("/Users/hbolfarine/Dropbox/DSS_MIX/DSS_MIX_06_2024/")
source("source/dcpossum_dens_comp.R")
source("source/dcpossum_sim_data_mix.R")
source("source/dcpossum_plots.R")
source("source/dcpossum_unc.R")
source("source/func_pred_laplace_temp.R")
source("source/dcpossum_clust.R")
# source("source/dcpossum_clust_compar.R")

# Numerical experiments 
# Acidity data 
# BP
#----------------------------------------------##----------------------------------------------#
#----------------------------------------------##----------------------------------------------#

BP.acidity = dcpossum.BP(acidity, kmax = 10, BP.run = F, quant.sample = 1000, data.bp = "acidity", pred.f = TRUE)

BP_comp_aci = plot.possum.uni(BP.acidity[[1]], kmax = 10, y.lim = c(-1,0.5))
BP_comp_aci

k.star = 4
possum.BP.acidity = possum.unc.quant.values(BP.acidity, K.sel = k.star)

plots.possum.exemp.BP = plot.possum.quant(possum.BP.acidity, K.sel = k.star, scale.plot = FALSE, 
                                       acidity, index.possum = FALSE, index.pred = "acidity", model = "BP")

possum_clust.BP = dc.possum.clust.uni(BP.acidity, acidity, K.sel = k.star, km = TRUE)

dat.BP = process_clustering_uni(BP.acidity[[6]], possum_clust.BP, k.star, acidity)

# y = -0.07, -0.15
plot.acidity.BP.3 = create_custom_plot(dat.BP, plots.possum.exemp.BP$dens.summ, k.star, y.c = -0.05, y.k = -0.1)

plot.acidity.BP.4 = create_custom_plot(dat.BP, plots.possum.exemp.BP$dens.summ, k.star, y.c = -0.05, y.k = -0.1)

post_comp_BP = plot_posterior_components(table(BP.acidity[[3]]))

BP_acidity = BP_comp_aci + plot.acidity.BP.3 + plot.acidity.BP.4 #+ plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)")))#+ post_comp_BP

p.acidity.BP.final <- ggdraw() +
  draw_plot(BP_acidity) +  # Your main plot
  draw_plot(plots.possum.exemp.BP$legend, x = 0.46, y = 0.74, width = 0.25, height = 0.25) +
  draw_plot(plots.possum.exemp.BP$legend, x = 0.78, y = 0.74, width = 0.25, height = 0.25) 

# ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/paper_test_overleaf/acidity_BP.pdf", 
#        plot = p.acidity.BP.final, width = 14.5, height = 4.5)

#----------------------------------------------##----------------------------------------------#
#----------------------------------------------##----------------------------------------------#

#DPM method
set.seed(1822)
# DPM.acidity = dcpossum.DPM(acidity, kmax = 10, quant.sample = 1000)
DPM.acidity = dcpossum.DPM.dir(acidity, kmax = 10, quant.sample = 1000, k0 = 1/10, pred.f = TRUE)

DPM_comp = plot.possum.uni(DPM.acidity[[1]], kmax = 10, y.lim = c(-1,0.5))
DPM_comp

k.star = 4
possum.DPM.acidity = possum.unc.quant.values(DPM.acidity, K.sel = k.star)

plots.possum.exemp.DPM = plot.possum.quant(possum.DPM.acidity, K.sel = k.star, acidity,
                                       index.possum = FALSE, index.pred = "acidity", model = "DPM")

possum_clust.DPM = dc.possum.clust.uni(DPM.acidity, acidity, K.sel = k.star, km = TRUE)

dat.DPM = process_clustering_uni(DPM.acidity[[6]], possum_clust.DPM, k.star, acidity)

# y = -0.07, -0.15
plot.acidity.DPM.3 = create_custom_plot(dat.DPM, plots.possum.exemp.DPM$dens.summ, k.star, y.c =  -0.05, y.k = -0.1)

plot.acidity.DPM.4 = create_custom_plot(dat.DPM, plots.possum.exemp.DPM$dens.summ, k.star, y.c =  -0.05, y.k = -0.1)

post_comp_DPM = plot_posterior_components(table(DPM.acidity[[3]]))

DPM_acidity = DPM_comp + plot.acidity.DPM.3 + plot.acidity.DPM.4 #+ plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)")))

p.acidity.DPM.final <- ggdraw() +
  draw_plot(DPM_acidity) +  # Your main plot
  draw_plot(plots.possum.exemp.DPM$legend ,x = 0.46, y = 0.74, width = 0.25, height = 0.25) +
  draw_plot(plots.possum.exemp.DPM$legend, x = 0.78, y = 0.74, width = 0.25, height = 0.25) 

# ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/paper_test_overleaf/acidity_DPM.pdf", 
#        plot = p.acidity.DPM.final, width = 14.5, height = 4.5)

#----------------------------------------------##----------------------------------------------#
#----------------------------------------------##----------------------------------------------#
# MFM method  
set.seed(1800) 
MFM.acidity = dcpossum.MFM(acidity, kmax = 10, quant.sample = 1000,pred.f = TRUE)

MFM_comp = plot.possum.uni(MFM.acidity[[1]], kmax = 10, y.lim = c(-1,0.6))
MFM_comp

k.star = 3
possum.MFM.acidity = possum.unc.quant.values(MFM.acidity, K.sel = k.star)

plots.possum.exemp.MFM = plot.possum.quant(possum.MFM.acidity, K.sel = k.star, scale.plot = FALSE, acidity, model = "MFM")

# MFM_comp + plots.possum.exemp.MFM$dens.summ

# # # Plot - number of components from the methods

possum_clust.MFM = dc.possum.clust.uni(MFM.acidity, acidity, K.sel = k.star, km = TRUE)

dat.MFM = process_clustering_uni(MFM.acidity[[6]], possum_clust.MFM, k.star, acidity)

# y = -0.07, -0.15
plot.acidity.MFM.3 = create_custom_plot(dat.MFM, plots.possum.exemp.MFM$dens.summ, k.star, y.c = -0.05, y.k = -0.1)

plot.acidity.MFM.4 = create_custom_plot(dat.MFM, plots.possum.exemp.MFM$dens.summ, k.star, y.c = -0.05, y.k = -0.1)

post_comp_MFM = plot_posterior_components(table(MFM.acidity[[3]]))

MFM_acidity = MFM_comp + plot.acidity.MFM.3 + plot.acidity.MFM.4 #+ plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)")))#+ post_comp_MFM

p.acidity.MFM.final <- ggdraw() +
  draw_plot(MFM_acidity) +  # Your main plot
  draw_plot(plots.possum.exemp.MFM$legend, x = 0.46, y = 0.74, width = 0.25, height = 0.25) +
  draw_plot(plots.possum.exemp.MFM$legend, x = 0.78, y = 0.74, width = 0.25, height = 0.25) 

# ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/acidity_MFM.pdf", 
#        plot = p.acidity.MFM.final, width = 14.5, height = 4.5)
# 
comp.all.acidity = plot.postcomp.uni(DPM.acidity,MFM.acidity,BP.acidity,K.true = FALSE)
# ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/acidity_post_all.pdf", 
#        plot = comp.all.acidity, width = 5, height = 3)

acidity.final.all = p.acidity.BP.final/p.acidity.DPM.final/p.acidity.MFM.final+ plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)")))

ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/sn-article-template/acidity_possum_all_post.pdf", 
       plot = acidity.final.all, width = 12.5, height = 11)


# comp.all.acidity + 

