setwd("/Users/hb23255/Dropbox/DSS_MIX/DSS_MIX_06_2024/")
source("source/dcpossum_dens_comp.R")
source("source/dcpossum_sim_data_mix.R")
source("source/dcpossum_plots.R")
source("source/dcpossum_unc.R")
source("source/func_pred_laplace_temp.R")
source("source/dcpossum_clust.R")
source("source/dcpossum_clust_compar.R")

set.seed(1820) 
y.data.1.1 = data_sim_func("example_02", n = 100)
# write.csv(y.data.1.1,"source/BP_source/sim_bp/bp_250/data_bp_01_01.csv")
y.data.1.2 = data_sim_func("example_02", n = 250)
# write.csv(y.data.1.2,"source/BP_source/sim_bp/bp_500/data_bp_01_02.csv")
y.data.1.3 = data_sim_func("example_02", n = 1000)
# write.csv(y.data.1.3,"source/BP_source/sim_bp/bp_1000/data_bp_01_03.csv")

# sim 
# BP.1 = dcpossum.BP(y.data.1.1, kmax = 10, BP.run = F, quant.sample = 1000, data.bp = "sim_250", pred.f = TRUE)
MFM.1 = dcpossum.MFM(y.data.1.1,kmax = 10, quant.sample = 1000, pred.f = TRUE)
DPM.1 = dcpossum.DPM(y.data.1.1,kmax = 10, quant.sample = 1000, pred.f = TRUE)

comp.1 = plot.possum.comp.uni(BP.1,DPM.1,MFM.1, title = "", kmax = 10, y.data = y.data.1.1, sel.K = 5)
comp.1

comp.2 = plot.postcomp.uni(DPM.1,MFM.1,BP.1,K.true = 5)
comp.2

possum.BP.1 = possum.unc.quant.values(BP.1, K.sel = 5)
possum.MFM.1 = possum.unc.quant.values(MFM.1, K.sel = 5)
possum.DPM.1 = possum.unc.quant.values(DPM.1, K.sel = 5)

plots.possum = plot.possum.unc.quant.dc(possum.BP.1,possum.MFM.1,possum.DPM.1, K.sel = 5, 
                                        y.data.1.1, index.possum = "example_02")
true.d = plots.possum$helling$td
# pred.d.BP.1 = BP.1[[8]]
# pred.d.MFM.1 = MFM.1[[8]]
# pred.d.DPM.1 = DPM.1[[8]]
#   
o.BP = plots.possum$helling$o.BP
m.BP = plots.possum$helling$m.BP
o.MFM = plots.possum$helling$o.MFM
m.MFM = plots.possum$helling$m.MFM
o.DPM = plots.possum$helling$o.DPM
m.DPM = plots.possum$helling$m.DPM

t.dens = true.d/sum(true.d)
# p.dens.BP = pred.d.BP.1/sum(pred.d.BP.1)
# p.dens.MFM = pred.d.MFM.1/sum(pred.d.MFM.1)
# p.dens.DPM = possum.DPM.1/sum(possum.DPM.1)

o.n.BP = o.BP/sum(o.BP)
o.m.BP = m.BP/sum(m.BP)
o.n.MFM = o.MFM/sum(o.MFM)
o.m.MFM = m.MFM/sum(m.MFM)
o.n.DPM = o.DPM/sum(o.DPM)
o.m.DPM = m.DPM/sum(m.DPM)

b1 = sqrt(sum((sqrt(t.dens) - sqrt(o.n.BP))^2)) / sqrt(2)
b2 = sqrt(sum((sqrt(t.dens) - sqrt(o.m.BP))^2)) / sqrt(2)

# b1p = sqrt(sum((sqrt(p.dens) - sqrt(o.n.BP))^2)) / sqrt(2)
# b2p = sqrt(sum((sqrt(p.dens) - sqrt(o.m.BP))^2)) / sqrt(2)

mfm1 = sqrt(sum((sqrt(t.dens) - sqrt(o.n.MFM))^2)) / sqrt(2)
mfm2 = sqrt(sum((sqrt(t.dens) - sqrt(o.m.MFM))^2)) / sqrt(2)

# mfm1p = sqrt(sum((sqrt(p.dens.MFM) - sqrt(o.n.MFM))^2)) / sqrt(2)
# mfm2p = sqrt(sum((sqrt(p.dens.MFM) - sqrt(o.m.MFM))^2)) / sqrt(2)

dpm1 = sqrt(sum((sqrt(t.dens) - sqrt(o.n.DPM))^2)) / sqrt(2)
dpm2 = sqrt(sum((sqrt(t.dens) - sqrt(o.m.DPM))^2)) / sqrt(2)

# dpm1p = sqrt(sum((sqrt(p.dens.DPM) - sqrt(o.n.DPM))^2)) / sqrt(2)
# dpm2p = sqrt(sum((sqrt(p.dens.DPM) - sqrt(o.m.DPM))^2)) / sqrt(2)

dat1 = data.frame(optimal  = c(b1,mfm1,dpm1), avdpossum = c(b2,mfm2,dpm2) , method = c("BP","MFM","DPM"), count.samp = as.factor(c(100,100,100)))

avg.possum.BP1 = BP.1[[1]]$avg.possum
avg.possum.MFM1 = MFM.1[[1]]$avg.possum
avg.possum.DPM1 = DPM.1[[1]]$avg.possum

plot(avg.possum.MFM1,type = "l")
lines(avg.possum.DPM1,type = "l", col = "blue")
lines(avg.possum.BP1,type = "l", col = "red")
abline(v = 5)

#-----------------------------##-----------------------------##-----------------------------##-----------------------------#
#-----------------------------##-----------------------------##-----------------------------##-----------------------------#

BP.12 = dcpossum.BP(y.data.1.2, kmax = 10, BP.run = F, quant.sample = 1000, data.bp = "sim_500", pred.f = TRUE)
MFM.12 = dcpossum.MFM(y.data.1.2,kmax = 10, quant.sample = 1000, pred.f = TRUE)
DPM.12 = dcpossum.DPM(y.data.1.2,kmax = 10, quant.sample = 1000, pred.f = TRUE)

comp.12 = plot.possum.comp.uni(BP.12,DPM.12,MFM.12, title = "", kmax = 10, y.data = y.data.1.2, sel.K = 5)
comp.12

comp.22 = plot.postcomp.uni(DPM.12,MFM.12,BP.12,K.true = 5)
comp.22

possum.BP.12 = possum.unc.quant.values(BP.12, K.sel = 5)
possum.MFM.12 = possum.unc.quant.values(MFM.12, K.sel = 5)
possum.DPM.12 = possum.unc.quant.values(DPM.12, K.sel = 5)

plots.possum2 = plot.possum.unc.quant.dc(possum.BP.12,possum.MFM.12,possum.DPM.12, K.sel = 5, 
                                        y.data.1.2, index.possum = "example_02")
true.d2 = plots.possum2$helling$td
o.BP2 = plots.possum2$helling$o.BP
m.BP2 = plots.possum2$helling$m.BP
o.MFM2 = plots.possum2$helling$o.MFM
m.MFM2 = plots.possum2$helling$m.MFM
o.DPM2 = plots.possum2$helling$o.DPM
m.DPM2 = plots.possum2$helling$m.DPM

t.dens2 = true.d2/sum(true.d2)

o.n.BP2 = o.BP2/sum(o.BP2)
o.m.BP2 = m.BP2/sum(m.BP2)
o.n.MFM2 = o.MFM2/sum(o.MFM2)
o.m.MFM2 = m.MFM2/sum(m.MFM2)
o.n.DPM2 = o.DPM2/sum(o.DPM2)
o.m.DPM2 = m.DPM2/sum(m.DPM2)

b12 = sqrt(sum((sqrt(t.dens2) - sqrt(o.n.BP2))^2)) / sqrt(2)
b22 = sqrt(sum((sqrt(t.dens2) - sqrt(o.m.BP2))^2)) / sqrt(2)

mfm12 = sqrt(sum((sqrt(t.dens2) - sqrt(o.n.MFM2))^2)) / sqrt(2)
mfm22 = sqrt(sum((sqrt(t.dens2) - sqrt(o.m.MFM2))^2)) / sqrt(2)

dpm12 = sqrt(sum((sqrt(t.dens2) - sqrt(o.n.DPM2))^2)) / sqrt(2)
dpm22 = sqrt(sum((sqrt(t.dens2) - sqrt(o.m.DPM2))^2)) / sqrt(2)

dat2 = data.frame(optimal  = c(b12,mfm12,dpm12), avdpossum = c(b22,mfm22,dpm22) , method = c("BP","MFM","DPM"), count.samp = as.factor(c(250,250,250)))

avg.possum.BP12 = BP.12[[1]]$avg.possum
avg.possum.MFM12 = MFM.12[[1]]$avg.possum
avg.possum.DPM12 = DPM.12[[1]]$avg.possum

plot(avg.possum.MFM12,type = "l")
lines(avg.possum.DPM12,type = "l", col = "blue")
lines(avg.possum.BP12,type = "l", col = "red")
abline(v = 5)

#-----------------------------##-----------------------------##-----------------------------##-----------------------------#
#-----------------------------##-----------------------------##-----------------------------##-----------------------------#

BP.13 = dcpossum.BP(y.data.1.3, kmax = 10, BP.run = F, quant.sample = 1000, data.bp = "sim_1000", pred.f = TRUE)
MFM.13 = dcpossum.MFM(y.data.1.3,kmax = 10, quant.sample = 1000, pred.f = TRUE)
DPM.13 = dcpossum.DPM(y.data.1.3,kmax = 10, quant.sample = 1000, pred.f = TRUE)

comp.13 = plot.possum.comp.uni(BP.13,DPM.13,MFM.13, title = "", kmax = 10, y.data = y.data.1.3, sel.K = 5)
comp.13

comp.23 = plot.postcomp.uni(DPM.13,MFM.13,BP.13,K.true = 5)
comp.23

possum.BP.13 = possum.unc.quant.values(BP.13, K.sel = 5)
possum.MFM.13 = possum.unc.quant.values(MFM.13, K.sel = 5)
possum.DPM.13 = possum.unc.quant.values(DPM.13, K.sel = 5)

plots.possum3 = plot.possum.unc.quant.dc(possum.BP.13,possum.MFM.13,possum.DPM.13, K.sel = 5, 
                                         y.data.1.3, index.possum = "example_02")

plots.possum3$dens.summ
plots.possum3$dens.summ
true.d3 = plots.possum3$helling$td
o.BP3 = plots.possum3$helling$o.BP
m.BP3 = plots.possum3$helling$m.BP
o.MFM3 = plots.possum3$helling$o.MFM
m.MFM3 = plots.possum3$helling$m.MFM
o.DPM3 = plots.possum3$helling$o.DPM
m.DPM3 = plots.possum3$helling$m.DPM

t.dens3 = true.d3/sum(true.d3)

o.n.BP3 = o.BP3/sum(o.BP3)
o.m.BP3 = m.BP3/sum(m.BP3)
o.n.MFM3 = o.MFM3/sum(o.MFM3)
o.m.MFM3 = m.MFM3/sum(m.MFM3)
o.n.DPM3 = o.DPM3/sum(o.DPM3)
o.m.DPM3 = m.DPM3/sum(m.DPM3)

b13 = sqrt(sum((sqrt(t.dens3) - sqrt(o.n.BP3))^2)) / sqrt(2)
b23 = sqrt(sum((sqrt(t.dens3) - sqrt(o.m.BP3))^2)) / sqrt(2)

mfm13 = sqrt(sum((sqrt(t.dens3) - sqrt(o.n.MFM3))^2)) / sqrt(2)
mfm23 = sqrt(sum((sqrt(t.dens3) - sqrt(o.m.MFM3))^2)) / sqrt(2)

dpm13 = sqrt(sum((sqrt(t.dens3) - sqrt(o.n.DPM3))^2)) / sqrt(2)
dpm23 = sqrt(sum((sqrt(t.dens3) - sqrt(o.m.DPM3))^2)) / sqrt(2)

dat3 = data.frame(optimal  = c(b13,mfm13,dpm13), avdpossum = c(b23,mfm23,dpm23) , method = c("BP","MFM","DPM"), count.samp = as.factor(c(1000,1000,1000)))

dat = data.frame(rbind(dat1, dat2, dat3))
dat.melt = melt(dat)
dat.melt$count.samp = as.numeric(as.character(dat.melt$count.samp ))
plot.hellinger = ggplot(dat.melt, aes(y = value, x = as.numeric(count.samp), color = method)) +
  geom_line() +
  facet_wrap(~variable) +
  labs(title="Hellinger Distance True density", x = "sample size")  


ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/Markdown_dssmix/anexos/hellinger.png", width = 10, height = 5)
  
avg.possum.BP13 = BP.13[[1]]$avg.possum
avg.possum.MFM13 = MFM.13[[1]]$avg.possum
avg.possum.DPM13 = DPM.13[[1]]$avg.possum

plot(avg.possum.MFM13,type = "l")
lines(avg.possum.DPM13,type = "l", col = "blue")
lines(avg.possum.BP13,type = "l", col = "red")
abline(v = 5)
 
#-----------------------------##-----------------------------##-----------------------------##-----------------------------#
#-----------------------------##-----------------------------##-----------------------------##-----------------------------#


# 
# possum.BP.1 = possum.unc.quant.values(BP.1, K.sel = 5)
# possum.MFM.1 = possum.unc.quant.values(MFM.1, K.sel = 5)
# possum.DPM.1 = possum.unc.quant.values(DPM.1, K.sel = 5)
# 
# 
# comp.1 = plot.possum.comp.uni(BP.1,DPM.1,MFM.1, title = "", kmax = 10, y.data = y.data.1, sel.K = 5)
# comp.1
# 
# # Plot - number of components from the methods
# comp.2 = plot.postcomp.uni(DPM.1,MFM.1,BP.1,K.true = 5)
# comp.2
# 
# #-----------------------------------------------------------------#
# # Load necessary libraries
# library(stats)
# 
# # Calculate the Hellinger distance
# # Load necessary libraries
# library(stats)
# 
# # Function to calculate Hellinger distance
# hellinger_distance <- function(p, q) {
#   if(length(p) != length(q)) {
#     stop("The lengths of the two distributions must be the same")
#   }
#   
#   # Calculate the Hellinger distance
#   distance <- sqrt(sum((sqrt(p) - sqrt(q))^2)) / sqrt(2)
#   return(distance)
# }
# 
# # Generate a sequence x on the reals
# set.seed(123)
# x <- rnorm(1000) # Example sequence from a normal distribution
# 
# # Estimate the empirical density from the sequence x
# empirical_density <- density(x)
# 
# # Define a theoretical density function, for example, normal distribution
# theoretical_density <- dnorm(empirical_density$x, mean=mean(x), sd=sd(x))
# 
# # Normalize the densities so that they sum to 1
# empirical_density$y <- empirical_density$y / sum(empirical_density$y)
# theoretical_density <- theoretical_density / sum(theoretical_density)
# 
# # Calculate the Hellinger distance
# distance <- hellinger_distance(empirical_density$y, theoretical_density)
# print(distance)
# 
