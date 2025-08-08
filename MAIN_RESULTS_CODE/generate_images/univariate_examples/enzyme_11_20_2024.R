setwd("/Users/hbolfarine/Dropbox/DSS_MIX/DSS_MIX_06_2024/")
source("source/dcpossum_dens_comp.R")
source("source/dcpossum_sim_data_mix.R")
source("source/dcpossum_plots.R")
source("source/dcpossum_unc.R")
source("source/func_pred_laplace_temp.R")
source("source/dcpossum_clust.R")
# source("source/dcpossum_clust_compar.R")

# Numerical experiments 
# Enzyme data 
load(file = "/Users/hbolfarine/Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/examples_data_BP/example_enzyme/enzyme.rda")
# write.csv(enzyme,"Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/generate_images/univariate_examples/enzime.csv", row.names = FALSE)
#----------------------------------------------##----------------------------------------------#
#----------------------------------------------##----------------------------------------------#
# BP
BP.enzyme = dcpossum.BP(enzyme, kmax = 10, BP.run = F, quant.sample = 1000, data.bp = "enzyme", pred.f = TRUE)

k.star = 2
BP_comp_sum.enzyme = plot.possum.uni(BP.enzyme[[1]], kmax = 10, y.lim = c(-1.8,0.35))
BP_comp_sum.enzyme

possum.BP.enzyme = possum.unc.quant.values(BP.enzyme, K.sel = k.star)

p.possum.exemp.enzyme.BP = plot.possum.quant(possum.BP.enzyme, K.sel = k.star, scale.plot = FALSE, 
                                          enzyme, index.possum = FALSE, index.pred = "enzyme", model = "BP")
p.possum.exemp.enzyme.BP$dens.summ

# BP_comp_sum.enzyme + p.possum.exemp.enzyme$dens.summ

possum_clust.BP = dc.possum.clust.uni(BP.enzyme, enzyme, K.sel = k.star, km = TRUE)

dat.BP = process_clustering_uni(BP.enzyme[[6]], possum_clust.BP, k.star, enzyme)

# y = -0.07, -0.15
plot.enzyme.BP.2 = create_custom_plot(dat.BP, p.possum.exemp.enzyme.BP$dens.summ, k.star, y.c = -0.12, y.k = -0.3)

plot.enzyme.BP.4 = create_custom_plot(dat.BP, p.possum.exemp.enzyme.BP$dens.summ, k.star, y.c = -0.12, y.k = -0.3)

post_comp_BP = plot_posterior_components(table(BP.enzyme[[3]]))

BP_enzyme_plot = BP_comp_sum.enzyme + plot.enzyme.BP.2 + plot.enzyme.BP.4 #+ plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)")))#+ post_comp_BP

p.enzyme.BP.final <- ggdraw() +
  draw_plot(BP_enzyme_plot) +  # Your main plot
  draw_plot(p.possum.exemp.enzyme.BP$legend, x = 0.36, y = 0.73, width = 0.25, height = 0.25) +
  draw_plot(p.possum.exemp.enzyme.BP$legend, x = 0.68, y = 0.73, width = 0.25, height = 0.25) 

# ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/enzyme_BP.pdf", 
#        plot = p.enzyme.BP.final, width = 14, height = 4.5)


#----------------------------------------------##----------------------------------------------#
#----------------------------------------------##----------------------------------------------#
set.seed(1822)
DPM.enzyme = dcpossum.DPM.dir(enzyme, kmax = 10, quant.sample = 1000, k0 = 1/10, pred.f = TRUE)

DPM_comp = plot.possum.uni(DPM.enzyme[[1]], kmax = 10, y.lim = c(-1.5,0.5))
DPM_comp

k.star = 2
possum.DPM.enzyme = possum.unc.quant.values(DPM.enzyme , K.sel = k.star)

plots.possum.exemp.DPM = plot.possum.quant(possum.DPM.enzyme , K.sel = k.star, scale.plot = FALSE, enzyme,
                                       index.possum = FALSE, model = "DPM")

# DPM_comp + plots.possum.exemp$dens.summ

possum_clust.DPM = dc.possum.clust.uni(DPM.enzyme, enzyme, K.sel = k.star, km = TRUE)

dat.DPM = process_clustering_uni(DPM.enzyme[[6]], possum_clust.DPM, k.star, enzyme)

# y = -0.07, -0.15
plot.enzyme.DPM.2 = create_custom_plot(dat.DPM, plots.possum.exemp.DPM$dens.summ, k.star, y.c = -0.12, y.k = -0.3)

plot.enzyme.DPM.4 = create_custom_plot(dat.DPM, plots.possum.exemp.DPM$dens.summ, k.star, y.c = -0.12, y.k = -0.3)

post_comp_DPM = plot_posterior_components(table(DPM.enzyme[[3]]))

DPM_enzyme_plot = DPM_comp + plot.enzyme.DPM.2 + plot.enzyme.DPM.4 #+ plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)")))#+ post_comp_DPM

p.enzyme.DPM.final <- ggdraw() +
  draw_plot(DPM_enzyme_plot) +  # Your main plot
  draw_plot(plots.possum.exemp.DPM$legend, x = 0.36, y = 0.73, width = 0.25, height = 0.25) +
  draw_plot(plots.possum.exemp.DPM$legend, x = 0.68, y = 0.73, width = 0.25, height = 0.25) 

# ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/enzyme_DPM.pdf", 
#        plot = p.enzyme.DPM.final, width = 14, height = 4.5)

#----------------------------------------------##----------------------------------------------#
#----------------------------------------------##----------------------------------------------#
set.seed(1822)
MFM.enzyme = dcpossum.MFM(enzyme, kmax = 10, quant.sample = 2000, pred.f = TRUE)

MFM_comp = plot.possum.uni(MFM.enzyme[[1]], kmax = 10)
MFM_comp

k.star = 4
possum.MFM.enzyme = possum.unc.quant.values(MFM.enzyme, K.sel = k.star)

plots.possum.exemp.MFM = plot.possum.quant(possum.MFM.enzyme, K.sel = k.star, scale.plot = FALSE, enzyme, model = "MFM")

# MFM_comp + plots.possum.exemp$dens.summ
plots.possum.exemp.MFM$dens.summ
# # # Plot - number of components from the methods

possum_clust.MFM = dc.possum.clust.uni(MFM.enzyme, enzyme, K.sel = k.star, km = TRUE)

dat.MFM = process_clustering_uni(MFM.enzyme[[6]], possum_clust.MFM, k.star, enzyme)

# y = -0.07, -0.15
plot.enzyme.MFM.2 = create_custom_plot(dat.MFM, plots.possum.exemp.MFM$dens.summ, k.star, y.c = -0.12, y.k = -0.3)

plot.enzyme.MFM.4 = create_custom_plot(dat.MFM, plots.possum.exemp.MFM$dens.summ, k.star, y.c = -0.12, y.k = -0.3)

post_comp_MFM = plot_posterior_components(table(MFM.enzyme[[3]]))

MFM_enzyme_plot = MFM_comp + plot.enzyme.MFM.2 + plot.enzyme.MFM.4 #+ plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)")))

p.enzyme.MFM.final <- ggdraw() +
  draw_plot(MFM_enzyme_plot) +  # Your main plot
  draw_plot(plots.possum.exemp.MFM$legend, x = 0.36, y = 0.73, width = 0.25, height = 0.25) +
  draw_plot(plots.possum.exemp.MFM$legend, x = 0.68, y = 0.73, width = 0.25, height = 0.25) 

# ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/paper_test_overleaf/enzyme_MFM.pdf", 
#        plot = p.enzyme.MFM.final, width = 14, height = 4.5)

# ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/enzyme_MFM_4.pdf", 
#        plot = MFM_enzyme_plot, width = 12, height = 3)

#----------------------------------------------##----------------------------------------------#
#----------------------------------------------##----------------------------------------------#

enzyme.final.all = p.enzyme.BP.final/p.enzyme.DPM.final/p.enzyme.MFM.final+ plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)")))

ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/sn-article-template/enzyme_possum_all_post.pdf", 
       plot = enzyme.final.all, width = 12.5, height = 11)


comp.all.enzyme = plot.postcomp.uni(DPM.enzyme,MFM.enzyme,BP.enzyme,K.true = FALSE)
ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/sn-article-template/enzyme_post_all.pdf", 
       plot = comp.all.enzyme, width = 5, height = 3)

comp.all.enzyme


