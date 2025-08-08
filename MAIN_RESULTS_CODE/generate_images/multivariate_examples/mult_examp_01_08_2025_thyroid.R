setwd("/Users/hbolfarine/Dropbox/DSS_MIX/DSS_MIX_06_2024/")
source("source/dcpossum_dens_comp.R")
source("source/dcpossum_sim_data_mix.R")
source("source/dcpossum_plots.R")
source("source/dcpossum_unc.R")
source("source/func_pred_laplace_temp.R")
source("source/dcpossum_clust.R")
# source("source/dcpossum_clust_compar.R")

#---------------------------------#
set.seed(12345)
y.data.thyroid = data_sim_func("thyroid")
y.data.clust = data_sim_func("thyroid_clust")

temp.MFM.mult = dcpossum.MFM.mult(y.data.thyroid, kmax = 10, quant.sample = 1000, 
                                  possum.samp = 1000)

temp.DPM.mult = dcpossum.DPM.dir.mult(y.data.thyroid, kmax = 10,
                                      possum.samp = 1000, f.pred = TRUE)

temp.SFM.mult = dcpossum.SFM.mult(y.data.thyroid, kmax = 10, col.data = 1:5, 
                                  clust.info = 0, quant.sample = 1000, post.size = 1000)

# Discrepancy plot
# temp.MFM.mult = temp.SFM.mult = temp.DPM.mult
comp.estim = plot.estim.comp.mult(temp.MFM.mult,temp.DPM.mult,temp.SFM.mult, title = "Discrenpancy function", sel.K = 3)
comp.estim

# Comparing posterior on the number of components
comp.model = plot.postcomp.mult(temp.MFM.mult,temp.DPM.mult,temp.SFM.mult, K.true = 3)
comp.model

#----------------------------------------------#

# comp.estim.DPM = plot.possum.uni(temp.DPM.mult[[1]], kmax = 10, sel.K = 3, y.lim = c(-12,2))
# comp.estim.DPM

# legend_title <- "Groups"
# 
# possum_clust = dc.possum.clust(temp.DPM.mult, y.data = y.data.thyroid, K.sel = 3, km = FALSE)
# optim.clust = dc.optim.clust(temp.DPM.mult, y.data.thyroid, K.sel = 3)
# 
# DPM.possum.clust = ggplot(possum_clust) +
#   geom_point(aes(x = dat.RT3U,y = dat.T4, col = clust.prob), size = 1.2) +
#   scale_color_distiller(palette = "Spectral", direction = -1, limits = c(0, 0.5)) +
#   theme_light() +
#   theme(legend.position="bottom") +
#   labs(title = "Possum uncertainty estimate", x = "RT3U", y = "T4")
# 
# DPM.optim.clust = ggplot(optim.clust) +
#   geom_point(aes(x = dat.RT3U,y = dat.T4, shape = as.factor(optim.clust), col = as.factor(optim.clust)), size = 1.2) +
#   theme_light() +
#   theme(legend.position="bottom") +
#   labs(title="Optimal Clustering classification", x = "RT3U", y = "T4") + scale_shape_manual(legend_title, values = c(15, 16, 17)) +
#   scale_colour_manual(legend_title, values=c("#619CFF", "#F8766D", "#00BA38"))
# 
# (comp.estim.DPM | DPM.optim.clust | DPM.possum.clust) + plot_annotation(tag_levels = list(c("a)", "b)", "c)", "d)")))
# ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/thyroid_DPM.pdf", width = 11, height = 4)

#-----------------------------------------------------#

# comp.estim.MFM = plot.possum.uni(temp.MFM.mult[[1]], kmax = 10, sel.K = 3, y.lim = c(-12,2))
# comp.estim.MFM

# legend_title <- "Groups"
# 
# possum_clust_MFM = dc.possum.clust(temp.MFM.mult, y.data = y.data.thyroid, K.sel = 3, km = FALSE)
# optim.clust_MFM = dc.optim.clust(temp.MFM.mult, y.data.thyroid, K.sel = 3)
# 
# MFM.possum.clust = ggplot(possum_clust_MFM) +
#   geom_point(aes(x = dat.RT3U,y = dat.T4, col = clust.prob), size = 1.2) +
#   scale_color_distiller(palette = "Spectral", direction = -1, limits = c(0, 0.5)) +
#   theme_light() +
#   theme(legend.position="bottom") +
#   labs(title = "Possum uncertainty estimate", x = "RT3U", y = "T4")
# 
# MFM.optim.clust = ggplot(optim.clust_MFM) +
#   geom_point(aes(x = dat.RT3U,y = dat.T4, shape = as.factor(optim.clust), col = as.factor(optim.clust)), size = 1.2) +
#   theme_light() +
#   theme(legend.position="bottom") +
#   labs(title="Optimal Clustering classification", x = "RT3U", y = "T4") + scale_shape_manual(legend_title, values = c(15, 16, 17)) +
#   scale_colour_manual(legend_title, values=c("#619CFF", "#F8766D", "#00BA38"))
# 
# (comp.estim.MFM | MFM.optim.clust | MFM.possum.clust) + plot_annotation(tag_levels = list(c("a)", "b)", "c)", "d)")))
# ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/thyroid_MFM.pdf", width = 11, height = 4)

#-----------------------------------------------------#

# comp.estim.SFM = plot.possum.uni(temp.SFM.mult[[1]], kmax = 10, sel.K = 3, y.lim = c(-12,2))
# comp.estim.SFM

# pairs(temp.DPM.mult[[6]])
# pairs(temp.DPM.mult[[4]][,,1])
# temp1 = Mclust(temp.DPM.mult[[6]],G = 3)
# temp2 = Mclust(y.pred,G = 3)
# # plot(temp1)
# temp1$parameters$variance$sigma[,,1]
# temp2$parameters$variance$sigma[,,2]
# temp1$parameters$mean
# tempclass1 = predict(temp1,newdata = y.data.thyroid)$classification
# tempclass2 = predict(temp2,newdata = y.data.thyroid)$classification
# 
# xtabs(~tempclass1 + tempclass2)


# temp.SFM.mult = temp.DPM.mult
# legend_title <- "Groups"
# possum_clust_SFM = possum_clust = data.clust
possum_clust_SFM = dc.possum.clust(temp.SFM.mult, y.data = y.data.thyroid, K.sel = 3, km = FALSE)
optim.clust_SFM = dc.optim.clust(temp.SFM.mult, y.data.thyroid, K.sel = 3)

# temp.DPM.mult
# legend_title <- "Groups"
# possum_clust_SFM = possum_clust = data.clust
possum_clust_DPM = dc.possum.clust(temp.DPM.mult, y.data = y.data.thyroid, K.sel = 3, km = FALSE)
optim.clust_DPM = dc.optim.clust(temp.DPM.mult, y.data.thyroid, K.sel = 3)


# pairs(possum_clust_DPM[,(1:5)],color = possum_clust_DPM$clust.prob)
# 
# library(GGally)
# ggpairs(possum_clust_DPM, aes(fill = clust.prob, alpha = 0.4))

clust.data.DPM = data.frame(possum_clust_DPM, clust.DPM = optim.clust_DPM$optim.clust)

DPM.possum.clust = ggplot(clust.data.DPM) +
  geom_point(aes(x = dat.RT3U,y = dat.T4, col = clust.prob, shape = as.factor(clust.DPM)), size = 1.5) +
  scale_color_distiller(palette = "Spectral", direction = -1, limits = c(0, 0.55)) +
  theme_light() +
  scale_shape_manual(name = "Groups", values=c(1,2,3)) + 
  theme(legend.position="bottom") +
  labs(title = "Possum uncertainty estimate", x = "RT3U", y = "T4")

clust.data.SFM = data.frame(possum_clust_SFM, clust.SFM = optim.clust_SFM$optim.clust)

SFM.possum.clust = ggplot(clust.data.SFM) +
  geom_point(aes(x = dat.RT3U,y = dat.T4, col = clust.prob, shape = as.factor(clust.SFM)), size = 1.5) +
  scale_color_distiller(palette = "Spectral", direction = -1, limits = c(0, 0.55)) +
  theme_light() +
  scale_shape_manual(name = "Groups", values=c(1,2,3)) + 
  # theme(legend.position="bottom") +
  theme(legend.position = "bottom", legend.title = element_text(size = 10), plot.title = element_text(size = 11, hjust = 0.5)) +
  labs(title = "Possum uncertainty estimate", x = "RT3U", y = "T4")

legend_title <- "Groups"
category_mapping <- c("Normal" = 2, "Hypo" = 1, "Hyper" = 3)
transformed_data <- thyroid %>%
  mutate(Category_Numeric = recode(Diagnosis, !!!category_mapping))

true_clust = ggplot(transformed_data) +
  geom_point(aes(x = RT3U, y = T4, shape = as.factor(Category_Numeric), col = as.factor(Category_Numeric)), size = 1.5) +
  theme_light() +
  theme(legend.position="bottom") +
  scale_shape_manual(name = "Groups", values=c(2,1,3)) + 
  labs(title="True", x ="RT3U", y = "T4") + #scale_shape_manual(legend_title, values = c(16, 17, 15)) +
  theme(legend.position = "bottom", legend.title = element_text(size = 10), plot.title = element_text(size = 11, hjust = 0.5)) +
  scale_colour_manual(legend_title, values=c("#619CFF", "#F8766D", "#00BA38"))

# theme(legend.position = "bottom", legend.title = element_text(size = 10), plot.title = element_text(size = 10, hjust = 0.5)) +

true_clust + SFM.possum.clust
# true_clust + DPM.possum.clust

clust.SFM = data.frame(transformed_data, SFM.clust = temp.SFM.mult[[9]])

clust_SFM = ggplot(clust.SFM) +
  geom_point(aes(x = RT3U, y = T4, shape = as.factor(SFM.clust), col = as.factor(SFM.clust)), size = 1.5) +
  theme_light() +
  theme(legend.position="bottom") +
  scale_shape_manual(name = "Groups", values=c(1,2,3,4,5,6)) + 
  labs(title="True", x ="RT3U", y = "T4") #+ #scale_shape_manual(legend_title, values = c(16, 17, 15)) +
  # scale_colour_manual(legend_title, values=c("#619CFF", "#F8766D", "#00BA38"))

true_clust + SFM.possum.clust + clust_SFM

# (comp.estim | comp.model) + plot_layout(widths = c(2, 1)) + plot_annotation(tag_levels = list(c("a)", "b)")))
p_comp_mult = ((comp.estim + comp.model) + plot_layout(widths = c(2, 1))) / (true_clust + SFM.possum.clust) + plot_annotation(tag_levels = list(c("(a)", "(b)","(c)","(d)")))
ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/thyroid_mult_all_pred.pdf", 
       plot = p_comp_mult, width = 11, height = 7)

#------------------------------------#

true_clust = ggplot(transformed_data) +
  geom_point(aes(x = RT3U, y = T4, shape = as.factor(Category_Numeric), col = as.factor(Category_Numeric)), size = 1.5) +
  theme_light() +
  theme(legend.position="bottom") +
  scale_shape_manual(name = "Groups", values=c(1,2,3)) + 
  labs(title="True", x ="RT3U", y = "T4") + #scale_shape_manual(legend_title, values = c(16, 17, 15)) +
  scale_colour_manual(legend_title, values=c("#619CFF", "#F8766D", "#00BA38"))

#------------------------------------#

SFM.possum.clust = ggplot(clust.data.SFM) +
  geom_point(aes(x = dat.RT3U,y = dat.T4, col = clust.prob, shape = as.factor(clust.SFM)), size = 2) +
  scale_color_distiller(palette = "Spectral", direction = -1, limits = c(0, max(clust.prob))) +
  theme_light() +
  scale_shape_manual(name = "Groups", values=c(1,2,3)) + 
  theme(legend.position="bottom") +
  labs(title = "Possum uncertainty estimate", x = "RT3U", y = "T4")

#------------------------------------#

SFM.optim.clust = ggplot(optim.clust_SFM) +
  geom_point(aes(x = dat.RT3U,y = dat.T4, shape = as.factor(optim.clust), col = as.factor(optim.clust)), size = 1.2) +
  theme_light() +
  theme(legend.position="bottom") +
  labs(title="Optimal Clustering classification", x = "RT3U", y = "T4") + scale_shape_manual(legend_title, values = c(16, 17, 15)) +
  scale_colour_manual(legend_title, values=c("#619CFF", "#F8766D", "#00BA38"))

category_mapping <- c("Normal" = 2, "Hypo" = 1, "Hyper" = 3)
transformed_data <- thyroid %>%
  mutate(Category_Numeric = recode(Diagnosis, !!!category_mapping))

true_clust = ggplot(transformed_data) +
  geom_point(aes(x = RT3U, y = T4, shape = as.factor(Category_Numeric), col = as.factor(Category_Numeric)), size = 1.2) +
  theme_light() +
  theme(legend.position="bottom") +
  labs(title="True", x ="RT3U", y = "T4") + scale_shape_manual(legend_title, values = c(16, 17, 15)) +
  scale_colour_manual(legend_title, values=c("#619CFF", "#F8766D", "#00BA38"))

# (15, 16, 17)
# (1, 2, 3)

# (2,3,1)
# (16, 17, 15)

table(optim.clust_SFM$optim.clust, transformed_data$Category_Numeric)
# table(optim.clust_SFM$optim.clust, thyroid$Diagnosis)


# (comp.estim | comp.model) / (true_clust + SFM.optim.clust + SFM.possum.clust) + plot_annotation(tag_levels = list(c("a)", "b)", "c)", "d)")))
(comp.estim | comp.model) + plot_layout(widths = c(2, 1)) + plot_annotation(tag_levels = list(c("a)", "b)")))
ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/thyroid_all_discr.pdf", width = 11, height = 3.5)

(true_clust + SFM.optim.clust + SFM.possum.clust) + plot_annotation(tag_levels = list(c("a)", "b)", "c)")))
ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/thyroid_SFM_true.pdf", width = 11, height = 4)


#-----------------------------------------------------#

# possum_clust.DPM = dc.possum.clust(temp.DPM.mult, y.data = y.data.thyroid, K.sel = 3, km = FALSE)
# 
# optim.clust.DPM = dc.optim.clust(temp.DPM.mult, y.data.thyroid, K.sel = 3)
# 
# ggplot(possum_clust.DPM) +
#   geom_point(aes(x = dat.RT3U,y = dat.T4, col = clust.prob), size = 1.2) +
#   scale_color_distiller(palette = "Spectral", direction = -1, limits = c(0, max(possum_clust$clust.prob))) +
#   theme_light() +
#   theme(legend.position="bottom") +
#   labs(title = "Possum uncertainty estimate", x = "RT3U", y = "T4")


# Insert table comparing classification

table(transformed_data$Category_Numeric, optim.clust.DPM$optim.clust)
classError(optim.clust$optim.clust, transformed_data$Category_Numeric)$errorRate
adjustedRandIndex(transformed_data$Category_Numeric, optim.clust.DPM$optim.clust)

ARI.SFM.orig = adjustedRandIndex(transformed_data$Category_Numeric, unlist(temp.SFM.mult[[9]]))
ARI.SFM.orig
classError(transformed_data$Category_Numeric, unlist(temp.SFM.mult[[9]]))$errorRate
# > ARI.SFM.orig
# [1] 0.8826543
# > classError(transformed_data$Category_Numeric, unlist(temp.SFM.mult[[9]]))$errorRate
# [1] 0.05581395

temp = Mclust(y.data.thyroid)
temp$classification

table(transformed_data$Category_Numeric, temp$classification)
classError(temp$classification, transformed_data$Category_Numeric)$errorRate
adjustedRandIndex(transformed_data$Category_Numeric, temp$classification)
# classError(temp$classification, transformed_data$Category_Numeric)$errorRate
# [1] 0.0372093
# > adjustedRandIndex(transformed_data$Category_Numeric, temp$classification)
# [1] 0.8771483


ARI.SFM.poss = adjustedRandIndex(transformed_data$Category_Numeric, optim.clust_SFM$optim.clust)
ARI.SFM.poss

classError(transformed_data$Category_Numeric, optim.clust_SFM$optim.clust)$errorRate


# > ARI.SFM.poss
# [1] 0.907901
# > classError(transformed_data$Category_Numeric, optim.clust_SFM$optim.clust)$errorRate
# [1] 0.02790698
