setwd("/Users/hbolfarine/Dropbox/DSS_MIX/DSS_MIX_06_2024/")
source("source/dcpossum_dens_comp.R")
source("source/dcpossum_sim_data_mix.R")
source("source/dcpossum_plots.R")
source("source/dcpossum_unc.R")
source("source/func_pred_laplace_temp.R")
source("source/dcpossum_clust.R")
# source("source/dcpossum_clust_compar.R")
 
set.seed(1822)
y.data.exam.03 = data_sim_func("example_mult_03_clust", n = 1000)

y.data.2 = y.data.exam.03[,c(1:2)]
y.clust = y.data.exam.03[,3]

temp.MFM.mult = dcpossum.MFM.mult(y.data.2, kmax = 10, 
                                  quant.sample = 1000, 
                                  possum.samp = 1000, f.pred = TRUE, pred.mod = "simulation")

# Working
temp.DPM.mult = dcpossum.DPM.dir.mult(y.data.2, kmax = 10,
                                  possum.samp = 1000, f.pred = TRUE, pred.mod = "simulation")

# Working
temp.SFM.mult = dcpossum.SFM.mult(y.data.2, kmax = 10, col.data = 1:2, 
                                  clust.info = 0, quant.sample = 1000, 
                                  possum.samp = 1000,
                                  post.size = 1000, f.pred = TRUE, pred.mod = "simulation")

# temp.MFM.mult = 
# temp.SFM.mult = temp.DPM.mult
# Number of components 
comp.estim = plot.estim.comp.mult(temp.MFM.mult,temp.DPM.mult,temp.SFM.mult, sel.K = 3, title = "Discrenpancy function")
comp.estim

# Number of components 
comp.mult = plot.postcomp.mult(temp.DPM.mult,temp.MFM.mult,temp.SFM.mult, K.true = 3)
comp.mult

MFM.dens = temp.MFM.mult
DPM.dens = temp.DPM.mult
SFM.dens = temp.SFM.mult

# Plot prediction 2d
plot.pred = plot_pred_dens_mult(SFM.dens,DPM.dens,MFM.dens, y.data.2, index = "example_02", sel.K = 3, plot_density = "simulation")
plot.pred

# Possum 
possum.SFM = possum.unc.quant.values.mult(temp.SFM.mult,K.sel = 3)
possum.DPM = possum.unc.quant.values.mult(temp.DPM.mult,K.sel = 3)
possum.MFM = possum.unc.quant.values.mult(temp.MFM.mult,K.sel = 3)

plot.possum.marg = plot.possum.uncertain.quant.dc.mult(possum.SFM, possum.MFM, possum.DPM, 
                                                       K.sel = 3, y.data.2, scale.plot = FALSE,
                                                       index.possum = FALSE, quant = c(0.025,0.975), 
                                                       examp = "simulation")

(comp.estim + comp.mult + plot_layout(widths = c(2, 1))) + plot_annotation(tag_levels = list(c("a)", "b)")))
ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/sim_plot_discr.pdf", width = 11, height = 3.5)

(plot.pred | plot.possum.marg) + plot_annotation(tag_levels = list(c("a)", "b)")))
ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/sim_plot_pred_post.pdf", width = 11, height = 4)

p_comp_pred = ((comp.estim + comp.mult) + plot_layout(widths = c(2, 1))) / (plot.pred + plot.possum.marg) + 
  plot_annotation(tag_levels = list(c("(a)", "(b)","(c)","(d)")))
ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/sim_mult_all_pred.pdf", 
                                                        plot = p_comp_pred, width = 11, height = 7.5) 
                                                 
# Plot possum cluster DPM
#--------------------------------------#

optim.clust = dc.optim.clust(temp.DPM.mult, y.data.2, K.sel = 3)
optim.clust.kmeans = dc.optim.clust.kmeans(temp.DPM.mult, y.data.2, K.sel = 3)
possum_clust = dc.possum.clust(temp.DPM.mult, y.data.2, K.sel = 3, km = TRUE)

clust_data = data.frame(possum_clust,optim.clust = as.factor(optim.clust$optim.clust), true.clust = as.factor(y.clust))

xtabs(~y.clust+optim.clust$optim.clust)

clust_data = clust_data %>% 
  mutate(rearange = case_when(
    optim.clust == 1 ~ 1,
    optim.clust == 2 ~ 2,
    optim.clust == 3 ~ 3
))

xtabs(~y.clust+clust_data$rearange)

plot.unc.clust = ggplot(clust_data ,aes(x = dat.V1, y = dat.V2, color = clust.prob, shape = as.factor(rearange))) +
  geom_point(size = 2) +
  theme_light() +
  xlab("") +
  ylab("") +
  scale_shape_manual(name = "Groups", values=c(1,2,3)) +
  theme(legend.position = "bottom", legend.title = element_text(size = 10), plot.title = element_text(size = 10, hjust = 0.5)) +
  # scale_color_distiller(palette = "Spectral", direction = -1, limits = c(0, .75)) 
  scale_color_distiller(name = "Probabilities", palette = "Spectral", direction = -1, limits = c(0, max(clust_data$clust.prob))) +
  labs(x ="y1", y = "y2") + #scale_shape_manual(legend_title, values = c(16, 17, 15)) +
  ggtitle("Probabilities Conditional summary")

# plot.true.clust = ggplot(clust_data) +
#   geom_point(aes(x = dat.V1, y = dat.V2, shape = as.factor(y.clust), col = as.factor(y.clust)), size = 2) +
#   theme_light() +
#   xlab("") +
#   ylab("") +
#   scale_shape_manual(name = "Clusters", values=c(1,2,3)) +
#   theme(legend.position = "bottom", legend.title = element_text(size = 10), plot.title = element_text(size = 10, hjust = 0.5)) +
#   labs(color = "Cluster Label", shape = "Cluster Label") +
#   ggtitle("True cluster labels")

legend_title = "Groups"
plot.true.clust = ggplot(clust_data) +
  geom_point(aes(x = dat.V1, y = dat.V2, shape = as.factor(y.clust), col = as.factor(y.clust)), size = 2) +
  theme_light() +
  theme(legend.position="bottom") +
  scale_shape_manual(name = "Groups", values=c(1,2,3)) +
  theme(legend.position = "bottom", legend.title = element_text(size = 10), plot.title = element_text(size = 10, hjust = 0.5)) +
  labs(title="True number of Groups", x ="y1", y = "y2") + #scale_shape_manual(legend_title, values = c(16, 17, 15)) +
  scale_colour_manual(legend_title, values=c("#619CFF", "#F8766D", "#00BA38"))

#----------------------------------------------------------------------#

plot.unc.km = ggplot(possum_clust,aes(x = dat.V1, y = dat.V2, color = clust.prob.km)) +
  geom_point(size = 0.5) +
  # geom_point() +
  theme_light() +
  xlab("") +
  ylab("") +
  theme(legend.position = "bottom", legend.title = element_text(size = 10), plot.title = element_text(size = 10, hjust = 0.5)) +
  # scale_color_distiller(palette = "Spectral", direction = -1, limits = c(0, .75)) 
  scale_color_distiller(name = "Probabilities", palette = "Spectral", direction = -1, limits = c(0, 0.75)) +
  ggtitle("Probabilities K-means summary")

plot.optim = ggplot(optim.clust,aes(x = dat.V1, y = dat.V2, shape = as.factor(optim.clust), color = as.factor(optim.clust))) +
  geom_point(size = 0.5) +
  # geom_point() +
  theme_light() +
  xlab("") +
  ylab("") +
  theme(legend.position = "bottom", legend.title = element_text(size = 10), plot.title = element_text(size = 10, hjust = 0.5)) +
  labs(color = "Cluster Label", shape = "Cluster Label") +
  ggtitle("Conditional summary cluster")

plot.optim.kmeans = ggplot(optim.clust.kmeans, aes(x = dat.V1, y = dat.V2, shape = as.factor(clust.km), color = as.factor(clust.km))) +
  geom_point(size = 0.5) +
  # geom_point() +
  theme_light() +
  xlab("") +
  ylab("") +
  labs(color = "Cluster Label", shape = "Cluster Label") +
  theme(legend.position = "bottom", legend.title = element_text(size = 10),plot.title = element_text(size = 10, hjust = 0.5)) +
  ggtitle("K-means summary cluster")

(plot.optim | plot.unc.clust | plot.optim.kmeans | plot.unc.km) + plot_annotation(tag_levels = list(c("a)", "b)", "c)", "d)")))

ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/sim_plot_DPM_sim.pdf", width = 11, height = 4)

plot.true.clust + plot.unc.clust + plot_annotation(tag_levels = list(c("(a)", "(b)")))
ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/sim_plot_DPM_sim.pdf", width = 10, height = 4.5)

# ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/sim_plot_DPM_sim_prob.pdf", 
#        plot = plot.unc.clust, width = 6, height = 4.5) 
# 
# ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/sim_plot_DPM_sim_true.pdf", 
#        plot = plot.true.clust, width = 6, height = 4.5) 


#--------------------------------------#
# Table
sim_Mclust = Mclust(y.data.2)

DPM.optim.clust = dc.optim.clust(temp.DPM.mult, y.data.2, K.sel = 3)
MFM.optim.clust = dc.optim.clust(temp.MFM.mult, y.data.2, K.sel = 3)
SFM.optim.clust = dc.optim.clust(temp.SFM.mult, y.data.2, K.sel = 3)

table(temp.MFM.mult[[9]][,1000])

clust.MFM = data.frame(clust.MFM = temp.MFM.mult[[9]][,1000])

clust.MFM = clust.MFM %>% 
  mutate(rearange = case_when(
    clust.MFM == 2 ~ 1,
    clust.MFM == 3 ~ 2,
    clust.MFM == 7 ~ 3
  ))

# table clustering
ARI.DPM.poss = adjustedRandIndex(y.clust, DPM.optim.clust$optim.clust)
ARI.DPM.poss
ARI.MFM.poss = adjustedRandIndex(y.clust, MFM.optim.clust$optim.clust)
ARI.MFM.poss
ARI.SFM.poss = adjustedRandIndex(y.clust, SFM.optim.clust$optim.clust)
ARI.SFM.poss
ARI.SFM.orig = adjustedRandIndex(y.clust, unlist(temp.SFM.mult[[9]]))
ARI.SFM.orig
ARI.DPM.orig = adjustedRandIndex(y.clust, temp.DPM.mult[[10]])
ARI.DPM.orig
ARI.MFM.orig = adjustedRandIndex(y.clust, clust.MFM$clust.MFM)
ARI.MFM.orig
ARI.Mclust = adjustedRandIndex(y.clust, sim_Mclust$classification)
ARI.Mclust

# > # table clustering
#   > ARI.DPM.poss = adjustedRandIndex(y.clust, DPM.optim.clust$optim.clust)
# > ARI.DPM.poss
# [1] 0.7713215
# > ARI.MFM.poss = adjustedRandIndex(y.clust, MFM.optim.clust$optim.clust)
# > ARI.MFM.poss
# [1] 0.7611891
# > ARI.SFM.poss = adjustedRandIndex(y.clust, SFM.optim.clust$optim.clust)
# > ARI.SFM.poss
# [1] 0.7708007
# > ARI.SFM.orig = adjustedRandIndex(y.clust, unlist(temp.SFM.mult[[9]]))
# > ARI.SFM.orig
# [1] 0.7622312
# > ARI.DPM.orig = adjustedRandIndex(y.clust, temp.DPM.mult[[10]])
# > ARI.DPM.orig
# [1] 0.5543276
# > ARI.MFM.orig = adjustedRandIndex(y.clust, clust.MFM$clust.MFM)
# > ARI.MFM.orig
# [1] 0.6489034
# > ARI.Mclust = adjustedRandIndex(y.clust, sim_Mclust$classification)
# > ARI.Mclust
# [1] 0.7600036


# category = numeric(1000)
# clust = numeric(1000)
# for(i in 1:1000){
#   if(length(unique(temp.MFM.mult[[9]][,i])) == 3){
#     category[i] = i
#   }else{
#     category[i] = 0
#   }
#   clust[i] = length(unique(temp.MFM.mult[[9]][,i]))
# }

classError(DPM.optim.clust$optim.clust, y.clust)$errorRate
classError(MFM.optim.clust$optim.clust, y.clust)$errorRate
classError(SFM.optim.clust$optim.clust, y.clust)$errorRate
classError(unlist(temp.SFM.mult[[9]]), y.clust)$errorRate
classError(clust.MFM$clust.MFM, y.clust)$errorRate
classError(temp.DPM.mult[[10]], y.clust)$errorRate
classError(sim_Mclust$classification, y.clust)$errorRate

# > classError(DPM.optim.clust$optim.clust, y.clust)$errorRate
# [1] 0.082
# > classError(MFM.optim.clust$optim.clust, y.clust)$errorRate
# [1] 0.086
# > classError(SFM.optim.clust$optim.clust, y.clust)$errorRate
# [1] 0.082
# > classError(unlist(temp.SFM.mult[[9]]), y.clust)$errorRate
# [1] 0.085
# > classError(clust.MFM$clust.MFM, y.clust)$errorRate
# [1] 0.133
# > classError(temp.DPM.mult[[10]], y.clust)$errorRate
# [1] 0.201
# > classError(sim_Mclust$classification, y.clust)$errorRate
# [1] 0.086

