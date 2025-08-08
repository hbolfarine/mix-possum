setwd("/Users/hb23255/Dropbox/DSS_MIX/DSS_MIX_06_2024/")
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

data_examp_02 = read.csv("/Users/hbolfarine/Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/examples_data_BP/example_02/data_bp_examp_02.csv")
data_examp_02 = data_examp_02$x

# examp_02_mclust = Mclust(data_examp_02)
# plot(examp_02_mclust)

BP.exmapl_02 = dcpossum.BP(data_examp_02, kmax = 10, BP.run = F, quant.sample = 1000, data.bp = "example_02", pred.f = TRUE)
# BP.exampl_02

# BP_comp_examp = plot.possum.uni(BP.exmapl_02[[1]], kmax = 10, y.data = data_examp_02, sel.K = 5)
BP_comp_examp = plot.possum.uni(BP.exmapl_02[[1]], kmax = 10, sel.K = 5, y.lim = c(-1,0.5))
BP_comp_examp

k.star = 5
possum.BP.example_02 = possum.unc.quant.values(BP.exmapl_02, K.sel = k.star)

plots.possum.exemp.BP = plot.possum.quant(possum.BP.example_02, K.sel = k.star, scale.plot = FALSE, 
                                       data_examp_02, index.possum = "example_02", index.pred = "example_02", model = "BP")

comp_BP = plot_posterior_components(table(BP.exmapl_02[[3]]), K.true = 5)

BP_examp = BP_comp_examp + plots.possum.exemp.BP$dens.summ + comp_BP #+ plot_annotation(tag_levels = list(c("a)", "b)", "c)")))

BP_examp.final <- ggdraw() +
  draw_plot(BP_examp) +  # Your main plot
  draw_plot(plots.possum.exemp.BP$legend, x = 0.465, y = 0.71, width = 0.25, height = 0.25)

ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/examp_possum_post_BP.pdf", 
       plot = BP_examp.final, width = 12.5, height = 3.5)


# possum_clust = dc.possum.clust.uni(BP.exmapl_02, data_examp_02, K.sel = k.star, km = TRUE)
# 
# dat = process_clustering_uni(BP.exmapl_02[[6]], possum_clust, k.star, data_examp_02)

# y = -0.07, -0.15
# plot.examp.BP = create_custom_plot(dat, plots.possum.exemp$dens.summ, k.star)

#----------------------------------------------##----------------------------------------------#
#----------------------------------------------##----------------------------------------------#

#DPM method
# data_examp_02 = data_sim_func("example_02", n = 2000)
data_examp_02 = read.csv("/Users/hbolfarine//Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/examples_data_BP/example_02/data_bp_examp_02.csv")
data_examp_02 = data_examp_02$x

DPM.examp_02 = dcpossum.DPM.dir(data_examp_02, kmax = 10, quant.sample = 1000, k0 = 1/5, pred.f = TRUE)

# DPM_comp = plot.possum.uni(DPM.examp_02[[1]], kmax = 10, y.data = data_examp_02, sel.K = 5)
DPM_comp = plot.possum.uni(DPM.examp_02[[1]], kmax = 10, sel.K = 5, y.lim = c(-1.2,0.5))
DPM_comp

k.star = 5
possum.DPM.example = possum.unc.quant.values(DPM.examp_02, K.sel = k.star)

plots.possum.exemp.DPM = plot.possum.quant(possum.DPM.example, K.sel = k.star, data_examp_02,
                                       index.possum = "example_02", index.pred = "example_02", model = "DPM")

comp_DPM = plot_posterior_components(table(DPM.examp_02[[3]]), K.true = 5)

DPM_examp = DPM_comp + plots.possum.exemp.DPM$dens.summ + comp_DPM  #+ plot_annotation(tag_levels = list(c("a)", "b)", "c)")))

DPM_examp.final <- ggdraw() +
  draw_plot(DPM_examp) +  # Your main plot
  draw_plot(plots.possum.exemp.DPM$legend, x = 0.465, y = 0.71, width = 0.25, height = 0.25)

# possum_clust = dc.possum.clust.uni(DPM.examp_02, data_examp_02, K.sel = k.star, km = TRUE)
# 
# dat = process_clustering_uni(DPM.examp_02[[6]], possum_clust, k.star, data_examp_02)
# 
# # y = -0.07, -0.15
# plot.example.DPM = create_custom_plot(dat, plots.possum.exemp$dens.summ, k.star)

ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/examp_possum_post_DPM.pdf", 
       plot = DPM_examp.final, width = 12.5, height = 3.5)


#----------------------------------------------##----------------------------------------------#
#----------------------------------------------##----------------------------------------------#
# MFM method  

MFM.examp_02 = dcpossum.MFM(data_examp_02, kmax = 10, quant.sample = 1000, pred.f = TRUE)

# MFM_comp = plot.possum.uni(MFM.examp_02[[1]], kmax = 10, y.data = data_examp_02, sel.K = 5)
MFM_comp = plot.possum.uni(MFM.examp_02[[1]], kmax = 10, sel.K = 5, y.lim = c(-1,0.5))
MFM_comp

k.star = 5
possum.MFM.example = possum.unc.quant.values(MFM.examp_02, K.sel = k.star)

plots.possum.exemp.MFM = plot.possum.quant(possum.MFM.example, K.sel = k.star, data_examp_02,
                                       index.possum = "example_02", index.pred = "example_02", model = "MFM")

comp_MFM = plot_posterior_components(table(MFM.examp_02[[3]]), K.true = 5)

MFM_examp = MFM_comp + plots.possum.exemp.MFM$dens.summ + comp_MFM #+ plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)")))

MFM_examp.final <- ggdraw() +
  draw_plot(MFM_examp) +  # Your main plot
  draw_plot(plots.possum.exemp.MFM$legend, x = 0.465, y = 0.71, width = 0.25, height = 0.25)

ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/examp_possum_post_MFM.pdf", 
       plot = MFM_examp.final, width = 12.5, height = 3.5)


examp.final.all = BP_examp.final/DPM_examp.final/MFM_examp.final + plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)"))) 

ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/sn-article-template/examp_possum_all_post_MFM.pdf", 
       plot = examp.final.all, width = 12.5, height = 11)

# possum_clust = dc.possum.clust.uni(MFM.examp_02, data_examp_02, K.sel = k.star, km = TRUE)
# 
# dat = process_clustering_uni(MFM.examp_02[[6]], possum_clust, k.star, data_examp_02)
# 
# # y = -0.07, -0.15
# plot.example.MFM = create_custom_plot(dat, plots.possum.exemp$dens.summ, k.star)
