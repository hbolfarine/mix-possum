setwd("/Users/hb23255/Dropbox/DSS_MIX/DSS_MIX_06_2024/")
source("source/dcpossum_dens_comp.R")
source("source/dcpossum_sim_data_mix.R")
source("source/dcpossum_plots.R")
source("source/dcpossum_unc.R")
source("source/func_pred_laplace_temp.R")
source("source/dcpossum_clust.R")
source("source/dcpossum_clust_compar.R")

set.seed(1822)
y.data.exam.03 = data_sim_func("example_mult_03_clust", n = 1000)
# plot(y.data.2[,1],y.data.2[,2])

y.data.2 = y.data.exam.03[,c(1:2)]
y.clust = y.data.exam.03[,3]

# Working
# Change 95% quantile
temp.MFM.mult = dcpossum.MFM.mult(y.data.2, kmax = 10, 
                                  quant.sample = 1000, 
                                  possum.samp = 1000, f.pred = TRUE)

y1.dist <- seq(-1, 12.5, length.out = 500)
y2.dist <- seq(-1, 12.5, length.out = 500)
# 
plot(y1.dist,temp.MFM.mult[[8]][,1])
plot(y2.dist,temp.MFM.mult[[8]][,2])

# Working
temp.DPM.mult = dcpossum.DPM.mult(y.data.2, kmax = 10, quant.sample = 1000,
                                  possum.samp = 1000, f.pred = TRUE)

plot(y1.dist,temp.DPM.mult[[8]][,1])
plot(y2.dist,temp.DPM.mult[[8]][,2])

# Working
temp.SFM.mult = dcpossum.SFM.mult(y.data.2, kmax = 10, col.data = 1:2, 
                                  clust.info = 0, quant.sample = 1000, 
                                  possum.samp = 1000,
                                  post.size = 1000, f.pred = TRUE)


plot(y1.dist,temp.SFM.mult[[8]][,1])
plot(y2.dist,temp.SFM.mult[[8]][,2])

# Change name of this plot - mult and uni
# change - tags in dcpossum - let insert tag order
# Plot discrepancy function given number of factors 

# temp.DPM.mult = temp.MFM.mult
# temp.MFM.mult = temp.DPM.mult
# temp.SFM.mult = temp.MFM.mult
comp.estim = plot.estim.comp.mult(temp.MFM.mult,temp.DPM.mult,temp.SFM.mult, sel.K = 3, title = "Discrenpancy function")
comp.estim

# change name of this plot
# change - tags in dcpossum - let insert tag order
# Plot number of factors 
comp.mult = plot.postcomp.mult(temp.DPM.mult,temp.MFM.mult,temp.SFM.mult, K.true = 3)
comp.mult

MFM.dens = temp.MFM.mult
DPM.dens = temp.DPM.mult
SFM.dens = temp.SFM.mult

plot.pred = plot_pred_dens_mult(SFM.dens,DPM.dens,MFM.dens, y.data.2, index = "example_02", sel.K = 3, plot_density = "simulation")
plot.pred

possum.SFM = possum.unc.quant.values.mult(temp.SFM.mult,K.sel = 3)
possum.DPM = possum.unc.quant.values.mult(temp.DPM.mult,K.sel = 3)
possum.MFM = possum.unc.quant.values.mult(temp.MFM.mult,K.sel = 3)

# add histogram 
# point estimate
plot.possum.marg = plot.possum.uncertain.quant.dc.mult(possum.SFM, possum.MFM, possum.DPM, 
                                                       K.sel = 3, y.data.2, scale.plot = FALSE,
                                                       index.possum = FALSE, quant = c(0.025,0.975),examp = "simulation")
plot.possum.marg

(comp.estim + comp.mult + plot_layout(widths = c(2, 1))) / (plot.pred | plot.possum.marg) + plot_annotation(tag_levels = list(c("a)", "b)", "c)", "d)")))
ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/sim_plot_06.pdf", width = 12, height = 8)

#---------------------------------#

possum_clust = dc.possum.clust(temp.SFM.mult, y.data.2, K.sel = 3, km = TRUE)

plot.unc.clust = ggplot(possum_clust,aes(x = dat.V1, y = dat.V2, color = clust.prob)) +
  geom_point(size = 0.6) +
  geom_point() +
  theme_light() +
  xlab("") +
  ylab("") +
  theme(legend.position="bottom") +
  # scale_color_distiller(palette = "Spectral", direction = -1, limits = c(0, .75)) 
  scale_color_distiller(palette = "Spectral", direction = -1, limits = c(0, 0.75)) +
  ggtitle("Uncertainty - soft")

plot.unc.km = ggplot(possum_clust,aes(x = dat.V1, y = dat.V2, color = clust.prob.km)) +
  geom_point(size = 0.6) +
  geom_point() +
  theme_light() +
  xlab("") +
  ylab("") +
  theme(legend.position="bottom") +
  # scale_color_distiller(palette = "Spectral", direction = -1, limits = c(0, .75)) 
  scale_color_distiller(palette = "Spectral", direction = -1, limits = c(0, 0.75)) +
  ggtitle("Uncertainty - K-means")

optim.clust = dc.optim.clust(temp.SFM.mult, y.data.2, K.sel = 3)

optim.clust.kmeans = dc.optim.clust.kmeans(temp.SFM.mult, y.data.2, K.sel = 3)

plot.optim = ggplot(optim.clust,aes(x = dat.V1, y = dat.V2, shape = as.factor(optim.clust), color = as.factor(optim.clust))) +
  geom_point(size = 0.6) +
  geom_point() +
  theme_light() +
  xlab("") +
  ylab("") +
  theme(legend.position="bottom") +
  ggtitle("Optimal - soft")

plot.optim.kmeans = ggplot(optim.clust.kmeans, aes(x = dat.V1, y = dat.V2, shape = as.factor(clust.km), color = as.factor(clust.km))) +
  geom_point(size = 0.6) +
  geom_point() +
  theme_light() +
  xlab("") +
  ylab("") +
  theme(legend.position="bottom") +
  ggtitle("Optimal - kmeans")

y.data.examp = as.data.frame(y.data.exam.03)

gb = GBclust.func(y.data.2, H = 3)

true_clust = ggplot(y.data.examp, aes(x = V1, y = V2, shape = as.factor(clust_assign), color = as.factor(clust_assign))) +
  geom_point(size = 0.6) +
  geom_point() +
  theme_light() +
  xlab("") +
  ylab("") +
  theme(legend.position="bottom") +
  ggtitle("True")

(true_clust + plot.optim + plot.optim.kmeans) / (plot.unc.clust + plot.unc.km + gb) 


(true_clust | plot.unc.clust | plot.unc.km | gb)
ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/sim_plot_09.pdf", width = 12, height = 4)

#---------------------------------#
set.seed(12345)
y.data.thyroid = data_sim_func("thyroid")
y.data.clust = data_sim_func("thyroid_clust")

temp.MFM.mult = dcpossum.MFM.mult(y.data.thyroid, kmax = 10, quant.sample = 1000, 
                                  possum.samp = 1000)

temp.DPM.mult = dcpossum.DPM.mult(y.data.thyroid, kmax = 10, quant.sample = 1000,
                                  possum.samp = 1000)

temp.SFM.mult = dcpossum.SFM.mult(y.data.thyroid, kmax = 10, col.data = 1:5, 
                                  clust.info = 0, quant.sample = 1000, post.size = 1000)


# temp.DPM.mult = temp.SFM.mult
# temp.MFM.mult = temp.SFM.mult
# temp.SFM.mult = temp.SFM.mult


comp.estim = plot.estim.comp.mult(temp.MFM.mult,temp.DPM.mult,temp.SFM.mult, title = "Discrenpancy function", sel.K = 3)
comp.estim

comp.model = plot.postcomp.mult(temp.MFM.mult,temp.DPM.mult,temp.SFM.mult, K.true = 3)
comp.model

possum_clust = dc.possum.clust(temp.SFM.mult, y.data = y.data.thyroid, K.sel = 3, km = FALSE)

optim.clust = dc.optim.clust(temp.SFM.mult, y.data.thyroid, K.sel = 3)
# 
# optim.clust.kmeans = dc.optim.clust.kmeans(temp.SFM.mult, y.data.thyroid, K.sel = 3)

category_mapping <- c("Normal" = 1, "Hypo" = 3, "Hyper" = 2)

# Transforming the categories using mutate and recode
transformed_data <- thyroid %>%
  mutate(Category_Numeric = recode(Diagnosis, !!!category_mapping))

clust.sfm = unlist(temp.SFM.mult[[9]])
kmax.sfm = max(clust.sfm)
legend_title <- "Groups"

colors.temp <- c("#619CFF", "#F8766D", "#00BA38", "#F564E3", "#00BFC4","#B79F00", "#F564E3", "#FF64B0", "#00BF7D", "#A3A500")
col.sfm = colors.temp[1:kmax.sfm]

shapes.temp <- c(15, 16, 17, 3, 7, 8, 18, 4, 19, 20)
shape.sfm = shapes.temp[1:kmax.sfm]

p1 = ggplot(transformed_data) +
  geom_point(aes(x = RT3U, y = T4, shape = as.factor(Category_Numeric), col = as.factor(Category_Numeric)), size = 1.2) +
  theme_light() +
  theme(legend.position="bottom") +
  labs(title="True", x ="RT3U", y = "T4") + scale_shape_manual(legend_title, values = c(15, 16, 17)) +
  scale_colour_manual(legend_title, values=c("#619CFF", "#F8766D", "#00BA38"))

optim.clust.temp = optim.clust
optim.clust.temp$optim.clust = as.factor(optim.clust.temp$optim.clust)
relabel_map <- c("2" = "Normal", "3" = "Hyper", "1" = "Hypo")
optim.clust.temp$optim.clust <- factor(sapply(optim.clust.temp$optim.clus, function(x) relabel_map[as.character(x)]))

# Make k-means

p2 = ggplot(optim.clust.temp) +
  geom_point(aes(x = dat.RT3U,y = dat.T4, shape = optim.clust, col = optim.clust), size = 1.2) +
  theme_light()  +
  theme(legend.position="bottom") +
  labs(title="Optimal Clustering classification -", x = "RT3U", y = "T4") + scale_shape_manual(legend_title, values = c(15, 16, 17)) +
  scale_colour_manual(legend_title, values=c("#619CFF", "#F8766D", "#00BA38"))

p3 = ggplot(possum_clust) +
  geom_point(aes(x = dat.RT3U,y = dat.T4, col = clust.prob), size = 1.2) +
  scale_color_distiller(palette = "Spectral", direction = -1, limits = c(0, max(possum_clust$clust.prob))) +
  theme_light() +
  theme(legend.position="bottom") +
  labs(title = "Possum uncertainty estimate", x = "RT3U", y = "T4")

p4 = ggplot(optim.clust) +
  geom_point(aes(x = dat.RT3U,y = dat.T4, shape = as.factor(clust.sfm), col = as.factor(clust.sfm)), size = 1.2) +
  theme_light()  +
  theme(legend.position="bottom") +
  labs(title="SFM classfication", x = "SRT3U", y = "T4") + scale_shape_manual(legend_title, values = shape.sfm) +
  scale_colour_manual(legend_title, values = col.sfm)

(comp.estim + comp.model + plot_layout(widths = c(2, 1))) / (p1 | p3 | p4)
ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/DSS_Mixture_models/sim_plot_10.pdf", width = 12, height = 8)


original_factor <- factor(c(1, 2, 3, 1, 2, 3, 1, 3, 2))

# Relabel the factors
relabel_map <- c("1" = "2", "2" = "3", "3" = "1")
new_factor <- factor(sapply(original_factor, function(x) relabel_map[as.character(x)]))

# Insert table comparing classification
