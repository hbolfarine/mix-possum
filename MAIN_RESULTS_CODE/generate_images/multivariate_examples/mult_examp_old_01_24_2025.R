setwd("/Users/hb23255/Dropbox/DSS_MIX/DSS_MIX_06_2024/")
source("source/dcpossum_dens_comp.R")
source("source/dcpossum_sim_data_mix.R")
source("source/dcpossum_plots.R")
source("source/dcpossum_unc.R")
source("source/func_pred_laplace_temp.R")
source("source/dcpossum_clust.R")
# source("source/dcpossum_clust_compar.R")
 
set.seed(1822)
y.data.old = data_sim_func("faithful")

# y.data.2 = y.data.exam.03[,c(1:2)]
# y.clust = y.data.exam.03[,3]

temp.MFM.mult = dcpossum.MFM.mult(y.data.old, kmax = 10, 
                                  quant.sample = 1000, 
                                  possum.samp = 1000, f.pred = TRUE, pred.mod = "old_faith")

# Working
temp.DPM.mult = dcpossum.DPM.dir.mult(y.data.old, kmax = 10,
                                      possum.samp = 1000, f.pred = TRUE, pred.mod = "old_faith")

# Working
temp.SFM.mult = dcpossum.SFM.mult(y.data.old, kmax = 10, col.data = 1:2, 
                                  clust.info = 0, quant.sample = 1000, 
                                  possum.samp = 1000,
                                  post.size = 1000, f.pred = TRUE, pred.mod = "old_faith")
# plot(temp.SFM.mult[[6]]$X1, temp.SFM.mult[[6]]$X2)

# temp = Mclust(temp.SFM.mult[[6]], G = 3, modelNames = "VVV")
# temp$parameters
# plot(temp)
# temp.MFM.mult = temp.SFM.mult = temp.DPM.mult
# Number of components 
comp.estim = plot.estim.comp.mult(temp.MFM.mult,temp.DPM.mult,temp.SFM.mult, title = "Discrepancy function", sel.K = 0, kmax = 10)
comp.estim

# Number of components 
comp.mult = plot.postcomp.mult(temp.DPM.mult,temp.MFM.mult,temp.SFM.mult, kmax = 10)
comp.mult

MFM.dens = temp.MFM.mult
DPM.dens = temp.DPM.mult
SFM.dens = temp.SFM.mult

# Plot prediction 2d
# y.data = y.data.old
plot.pred = plot_pred_dens_mult(SFM.dens,DPM.dens,MFM.dens, y.data.old, index = 1, sel.K = 2, plot_density = "applic")
plot.pred

# Possum 
possum.SFM = possum.unc.quant.values.mult(temp.SFM.mult,K.sel = 2)
possum.DPM = possum.unc.quant.values.mult(temp.DPM.mult,K.sel = 2)
possum.MFM = possum.unc.quant.values.mult(temp.MFM.mult,K.sel = 2)

plot.possum.marg = plot.possum.uncertain.quant.dc.mult(possum.SFM, possum.MFM, possum.DPM, 
                                                       K.sel = 2, y.data.2, scale.plot = FALSE,
                                                       index.possum = FALSE, quant = c(0.025,0.975), 
                                                       examp = "old_faith")

p_comp_old = (comp.estim + comp.mult + plot_layout(widths = c(2, 1))) + plot_annotation(tag_levels = list(c("a)", "b)")))
# ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/sim_plot_discr.pdf", width = 11, height = 3.5)
# 
ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/old_mult_all.pdf", 
       plot = p_comp_old, width = 11, height = 3.5)


p_pred_old = (plot.pred + plot.possum.marg) + plot_layout(widths = c(1.7, 2)) + plot_annotation(tag_levels = list(c("a)", "b)")))
# ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/sim_plot_pred_post.pdf", width = 11, height = 4)
ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/old_plot_pred_post.pdf", 
       plot = p_pred_old, width = 11, height = 4.5)

# Plot possum cluster DPM
#--------------------------------------#

p_comp_pred = ((comp.estim + comp.mult) + plot_layout(widths = c(2, 1))) / (plot.pred + plot.possum.marg) + plot_annotation(tag_levels = list(c("a)", "b)","c)","d)")))
ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/old_mult_all_pred.pdf", 
       plot = p_comp_pred, width = 11, height = 8)


optim.clust = dc.optim.clust(temp.DPM.mult, y.data.old, K.sel = 2)
optim.clust.kmeans = dc.optim.clust.kmeans(temp.SFM.mult, y.data.old, K.sel = 2)
possum_clust = dc.possum.clust(temp.DPM.mult, y.data.old, K.sel = 2, km = TRUE)

plot.unc.clust = ggplot(possum_clust,aes(x = y.data.old[,1], y = y.data.old[,2], color = clust.prob, shape = as.factor(optim.clust$optim.clust))) +
  geom_point(size = 0.5) +
  geom_point() +
  theme_light() +
  xlab("Eruptions") +
  ylab("Waiting") +
  scale_shape_manual(name = "Clusters", values=c(3,4)) +
  theme(legend.position = "bottom", legend.title = element_text(size = 10), plot.title = element_text(size = 10, hjust = 0.5)) +
  # scale_color_distiller(palette = "Spectral", direction = -1, limits = c(0, .75)) 
  scale_color_distiller(name = "Probabilities", palette = "Spectral", direction = -1, limits = c(0, 0.5)) +
  ggtitle("Probabilities Conditional summary")

plot.unc.km = ggplot(possum_clust,aes(x = y.data.old[,1], y = y.data.old[,2],shape = as.factor(optim.clust.kmeans$clust.km), color = clust.prob.km)) +
  geom_point(size = 0.5) +
  geom_point() +
  theme_light() +
  xlab("Eruptions") +
  ylab("Waiting") +
  scale_shape_manual(name = "Clusters", values=c(3,8,4)) +
  theme(legend.position = "bottom", legend.title = element_text(size = 10), plot.title = element_text(size = 10, hjust = 0.5)) +
  # scale_color_distiller(palette = "Spectral", direction = -1, limits = c(0, .75)) 
  scale_color_distiller(name = "Probabilities", palette = "Spectral", direction = -1, limits = c(0, 0.5)) +
  ggtitle("Probabilities k-means summary")

p_clust_old = plot.unc.clust + plot.unc.km + plot_annotation(tag_levels = list(c("a)", "b)")))
ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/old_plot_pred_clust.pdf", 
       plot = p_clust_old, width = 11, height = 4.5)


library(mclust)
temp = Mclust(y.data.old,modelNames = "VVV",G = 3)
plot(temp)

optim.clust %>% 
  arrange(desc(dat.waiting))

# plot.optim = ggplot(optim.clust,aes(x = y.data.old[,1], y = y.data.old[,2], shape = as.factor(optim.clust), color = as.factor(optim.clust))) +
#   geom_point(size = 0.5) +
#   geom_point() +
#   theme_light() +
#   xlab("") +
#   ylab("") +
#   theme(legend.position = "bottom", legend.title = element_text(size = 10), plot.title = element_text(size = 10, hjust = 0.5)) +
#   labs(color = "Cluster Label", shape = "Cluster Label") +
#   ggtitle("Conditional summary cluster")
# 
# plot.optim.kmeans = ggplot(optim.clust.kmeans, aes(x = y.data.old[,1], y = y.data.old[,2], shape = as.factor(clust.km), color = as.factor(clust.km))) +
#   geom_point(size = 0.5) +
#   geom_point() +
#   theme_light() +
#   xlab("") +
#   ylab("") +
#   labs(color = "Cluster Label", shape = "Cluster Label") +
#   theme(legend.position = "bottom", legend.title = element_text(size = 10),plot.title = element_text(size = 10, hjust = 0.5)) +
#   ggtitle("K-means summary cluster")

# (plot.unc.clust | plot.unc.km) + plot_annotation(tag_levels = list(c("a)", "b)", "c)", "d)")))

# ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/sim_plot_DPM_sim.pdf", width = 11, height = 4)

#--------------------------------------#
# Table
# sim_Mclust = Mclust(y.data.old)
# 
# DPM.optim.clust = dc.optim.clust(temp.DPM.mult, y.data.old, K.sel = 3)
# MFM.optim.clust = dc.optim.clust(temp.MFM.mult, y.data.old, K.sel = 3)
# SFM.optim.clust = dc.optim.clust(temp.SFM.mult, y.data.old, K.sel = 3)
# 
# # table clustering
# ARI.DPM.poss = adjustedRandIndex(y.clust, DPM.optim.clust$optim.clust)
# ARI.MFM.poss = adjustedRandIndex(y.clust, MFM.optim.clust$optim.clust)
# ARI.SFM.poss = adjustedRandIndex(y.clust, SFM.optim.clust$optim.clust)
# ARI.SFM.orig = adjustedRandIndex(y.clust, unlist(temp.SFM.mult[[9]]))
# ARI.DPM.orig = adjustedRandIndex(y.clust, temp.DPM.mult[[9]][,998])
# ARI.MFM.orig = adjustedRandIndex(y.clust, temp.MFM.mult[[9]][,982])
# ARI.Mclust = adjustedRandIndex(y.clust, sim_Mclust$classification)
# 
# # category = numeric(1000)
# # clust = numeric(1000)
# # for(i in 1:1000){
# #   if(length(unique(temp.MFM.mult[[9]][,i])) == 3){
# #     category[i] = i
# #   }else{
# #     category[i] = 0
# #   }
# #   clust[i] = length(unique(temp.MFM.mult[[9]][,i]))
# # }
# 
# classError(DPM.optim.clust$optim.clust, y.clust)$errorRate
# classError(MFM.optim.clust$optim.clust, y.clust)$errorRate
# classError(SFM.optim.clust$optim.clust, y.clust)$errorRate
# classError(unlist(temp.SFM.mult[[9]]), y.clust)$errorRate
# classError(temp.DPM.mult[[9]][,998], y.clust)$errorRate
# classError(temp.DPM.mult[[9]][,982], y.clust)$errorRate
# classError(sim_Mclust$classification, y.clust)$errorRate


