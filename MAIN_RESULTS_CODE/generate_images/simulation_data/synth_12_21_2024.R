setwd("/Users/hb23255/Dropbox/DSS_MIX/DSS_MIX_06_2024/")
source("source/dcpossum_dens_comp.R")
source("source/dcpossum_sim_data_mix.R")
source("source/dcpossum_plots.R")
source("source/dcpossum_unc.R")
source("source/func_pred_laplace_temp.R")
source("source/dcpossum_clust.R")
source("source/dcpossum_clust_compar.R")


# Numerical experiments 

# Univariate models - desnity estimation

# Run with synthetic data

setwd("/Users/hb23255/Dropbox/DSS_MIX/DSS_MIX_06_2024/")
y.data.1 = read.csv("source/BP_source/examples_data_BP/example_02/data_bp_examp_02.csv")
y.data.1 = y.data.1$x

#----------------------------------------------#

set.seed(1820) 
y.data.examp = data_sim_func("example_02", n = 1000)
y.data.examp = y.data.1

DPM.examp = dcpossum.DPM.dir(y.data.examp, kmax = 10, quant.sample = 1000, k0 = 0.1, pred.f = TRUE)

K.star = 5

DPM_comp = plot.possum.uni(DPM.examp[[1]], kmax = 10, sel.K = K.star, y.lim = c(-0.9,0.4))
DPM_comp

possum.DPM.examp = possum.unc.quant.values(DPM.examp, K.sel = K.star)

plots.possum.exemp.DPM = plot.possum.quant(possum.DPM.examp, K.sel = K.star, y.data.examp,
                                       index.possum = "example_02", index.pred = "example_02")

DPM_comp + plots.possum.exemp.DPM$dens.summ

# possum_clust = dc.possum.clust.uni(DPM.examp, y.data.examp, K.sel = K.star, km = TRUE)
# 
# dat = process_clustering_uni(DPM.examp[[6]], possum_clust, K.star, y.data.examp)
# 
# # y = -0.07, -0.15
# plot.examp.DPM = create_custom_plot(dat, plots.possum.exemp$dens.summ, K.star)
# 
# table(DPM.examp[[3]])

DPM_post_clust = plot_posterior_components(table(DPM.examp[[3]]))

DPM_model_eval = DPM_comp + plots.possum.exemp.DPM$dens.summ + DPM_post_clust

ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/examp_possum_post_DPM.pdf", 
       plot = DPM_model_eval, width = 15, height = 3.5)

#----------------------------------------------#

set.seed(1820) 
# y.data.examp = data_sim_func("example_02", n = 1000)

MFM.examp = dcpossum.MFM(y.data.examp, kmax = 10, quant.sample = 5000)

K.star = 5

MFM_comp = plot.possum.uni(MFM.examp[[1]], kmax = 10, y.data = y.data.examp, sel.K = K.star,  y.lim = c(-0.9,0.4))
MFM_comp

possum.MFM.examp = possum.unc.quant.values(MFM.examp, K.sel = K.star)

plots.possum.exemp.MFM = plot.possum.quant(possum.MFM.examp, K.sel = K.star, y.data.examp,
                                       index.possum = "example_02", index.pred = "example_02")

MFM_comp + plots.possum.exemp.MFM$dens.summ

# possum_clust = dc.possum.clust.uni(DPM.examp, y.data.examp, K.sel = K.star, km = TRUE)
# 
# dat = process_clustering_uni(DPM.examp[[6]], possum_clust, K.star, y.data.examp)

# y = -0.07, -0.15
# plot.examp.DPM = create_custom_plot(dat, plots.possum.exemp$dens.summ, K.star)

# table(DPM.examp[[3]])

MFM_post_clust = plot_posterior_components(table(MFM.examp[[3]]))

MFM_model_eval = MFM_comp + plots.possum.exemp.MFM$dens.summ + MFM_post_clust

ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/examp_possum_post_MFM.pdf", 
       plot = MFM_model_eval, width = 15, height = 3.5)

#----------------------------------------------#

BP.mix = dcpossum.BP(y.data.1, kmax = 10, BP.run = F, quant.sample = 5000, data.bp = "example_02", pred.f = FALSE, scale.plot = c(10,40))

# hist(BP.mix[[6]],breaks = 100)

BP_comp_sum = plot.possum.uni(BP.mix[[1]], kmax = 10, y.data = y.data.1, sel.K = 5,  y.lim = c(-0.9,0.4))
BP_comp_sum

possum.BP.1 = possum.unc.quant.values(BP.mix, K.sel = 5)

plots.possum.exemp.BP = plot.possum.quant(possum.BP.1, K.sel = 5, y.data.1,
                                       index.possum = "example_02", index.pred = "example_02")
plots.possum.exemp.BP$dens.summ

BP_comp_sum + plots.possum.exemp.BP$dens.summ

BP_post_clust = plot_posterior_components(table(BP.mix[[3]]))

BP_model_eval = BP_comp_sum + plots.possum.exemp.BP$dens.summ + BP_post_clust

ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/examp_possum_post_BP.pdf", 
       plot = BP_model_eval, width = 15, height = 3.5)

