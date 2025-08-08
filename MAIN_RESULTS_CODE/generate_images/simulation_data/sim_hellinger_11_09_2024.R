setwd("/Users/hbolfarine/Dropbox/DSS_MIX/DSS_MIX_06_2024/")
source("source/dcpossum_dens_comp.R")
source("source/dcpossum_sim_data_mix.R")
source("source/dcpossum_plots.R")
source("source/dcpossum_unc.R")
source("source/func_pred_laplace_temp.R")
source("source/dcpossum_clust.R")
# source("source/dcpossum_clust_compar.R")

set.seed(1820) 
y.data.1.1 = data_sim_func("example_02", n = 100)
# write.csv(y.data.1.1,"source/BP_source/sim_bp/bp_250/data_bp_01_01.csv")
y.data.1.2 = data_sim_func("example_02", n = 250)
# write.csv(y.data.1.2,"source/BP_source/sim_bp/bp_500/data_bp_01_02.csv")
y.data.1.3 = data_sim_func("example_02", n = 1000)
# write.csv(y.data.1.3,"source/BP_source/sim_bp/bp_1000/data_bp_01_03.csv")

# y.data.1.1 = y.data.1.3

# sim 
# BP.1 = dcpossum.BP(y.data.1.1, kmax = 10, BP.run = F, quant.sample = 1000, data.bp = "sim_250", pred.f = TRUE)
# MFM.1 = dcpossum.MFM(y.data.1.1,kmax = 10, quant.sample = 1000, pred.f = TRUE)
DPM.1 = dcpossum.DPM.dir(data_examp_02, kmax = 10, quant.sample = 1000, k0 = 1/5, pred.f = TRUE)

# comp.1 = plot.possum.comp.uni(BP.1,DPM.1,MFM.1, title = "", kmax = 10, y.data = y.data.1.1, sel.K = 5)
# comp.1

# comp.2 = plot.postcomp.uni(DPM.1,MFM.1,BP.1,K.true = 5)
# comp.2
MFM_comp_sum = plot.possum.uni(MFM.1[[1]], kmax = 10, y.data = y.data.1.1, sel.K = 5)
MFM_comp_sum

DPM_comp_sum = plot.possum.uni(DPM.1[[1]], kmax = 10, y.data = y.data.1.1, sel.K = 5)
DPM_comp_sum

possum.MFM.1 = possum.unc.quant.values(MFM.1, K.sel = 5)
possum.DPM.1 = possum.unc.quant.values(DPM.1, K.sel = 5)

# plots.possum = plot.possum.unc.quant.dc(possum.BP.1,possum.MFM.1,possum.DPM.1, K.sel = 5, 
#                                         y.data.1.1, index.possum = "example_02")

plots.possum.MFM = plot.possum.quant(possum.MFM.1, K.sel = 5, scale.plot = FALSE, 
                                     y.data.1.1, index.possum = "example_02")

plots.possum.DPM = plot.possum.quant(possum.DPM.1, K.sel = 5, scale.plot = FALSE, 
                                     y.data.1.1, index.possum = "example_02")
plots.possum.MFM$dens.summ
plots.possum.DPM$dens.summ

# library(mclust)
# model <- Mclust(y.data.1.1)
# summary(model)
# 
# bic_value <- model$bic
# print(bic_value)
# 
# plot(model)

true.d = plots.possum.MFM$helling$td
t.dens = true.d/sum(true.d)

optim = plots.possum.MFM$helling$optim.pred
optim.pred = optim/sum(optim)

mean.pred = plots.possum.MFM$helling$mean.pred
mean.possum = mean.pred/sum(mean.pred)

h1_optim = sqrt(sum((sqrt(t.dens) - sqrt(optim.pred))^2)) / sqrt(2)
h2_mean = sqrt(sum((sqrt(t.dens) - sqrt(mean.possum))^2)) / sqrt(2)


dat1 = data.frame(optimal  = c(b1,mfm1,dpm1), avdpossum = c(b2,mfm2,dpm2) , method = c("BP","MFM","DPM"), count.samp = as.factor(c(100,100,100)))

plot(avg.possum.MFM1,type = "l")
lines(avg.possum.DPM1,type = "l", col = "blue")
lines(avg.possum.BP1,type = "l", col = "red")
abline(v = 5)

# #-----------------------------##-----------------------------##-----------------------------##-----------------------------#
# #-----------------------------##-----------------------------##-----------------------------##-----------------------------#

results <- data.frame(
  method = character(),
  sample_size = integer(),
  iteration = integer(),
  h1_optim = numeric(),
  h2_mean = numeric(),
  exepec_post = numeric(),
  stringsAsFactors = FALSE
)

sample.sim = c(100,250,1000)
# method = c("MFM","DPM")
method = c("DPM")
iterations = 50
for(met in 1:length(method)){
  
  for(samp in 1:length(sample.sim)){
    
    for(iter in 1:iterations){
      
      y.data = data_sim_func("example_02", n = sample.sim[samp])
      
      # func_name = paste("dcpossum", method[j], sep = ".")
      # 
      # out.optim = do.call(func_name, list(y.data, kmax = 10, quant.sample = 1000, pred.f = TRUE))
      
      out.optim = dcpossum.DPM.dir(y.data, kmax = 10, quant.sample = 1000, k0 = 1/5, pred.f = TRUE)
      # comp_sum = plot.possum.uni(out.optim[[1]], kmax = 10, y.data = y.data, sel.K = 5)
      
      possum = possum.unc.quant.values(out.optim, K.sel = 5)
      
      # possum_quant = plot.possum.quant(possum, K.sel = 5, scale.plot = FALSE, 
      #                                  y.data, index.possum = "example_02")
      
      #-----------------------------------------------#
      
      # True
      true.d =  0.2*dnorm(y.data,mean = 19,sd = sqrt(5)) + 0.2*dnorm(y.data,mean = 19,sd = 1) +
        0.25*dnorm(y.data,mean = 23, sd = 1) + 0.2*dnorm(y.data,mean = 29, sd = 1) + 0.15*dnorm(y.data,mean = 33, sd = sqrt(2))
      
      # Pred point
      predi.point = numeric(length(y.data))
      prob.point = possum$point.estim$pro[[5]]
      mean.point = possum$point.estim$mean[[5]]
      s2.point = possum$point.estim$sigmasq[[5]]
      for(i in 1:length(y.data)){
        predi.point[i] = sum(prob.point*dnorm(y.data[i],mean = mean.point,sd = sqrt(s2.point)))
      }
      
      # Possum Expec
      list.param.method = list(mean = possum[[3]]$mean_post, 
                               sigmasq = possum[[3]]$sigmasq_post, 
                               pro = possum[[3]]$prob_post)
      
      N = nrow(list.param.method$mean)
      
      predi.mat = matrix(0,N,length(y.data))
      for(i in 1:N){
        for(j in 1:length(y.data)){
          predi.mat[i,j] = sum(list.param.method$pro[i,]*dnorm(y.data[j],mean = list.param.method$mean[i,],sd = sqrt(list.param.method$sigmasq[i,])))
        }
      }
      
      mean.pred.possum = apply(predi.mat, 2, mean)
      
      #-----------------------------------------------#
      
      t.dens = true.d/sum(true.d)
      
      optim.pred = predi.point/sum(predi.point)
      
      mean.possum = mean.pred.possum/sum(mean.pred.possum)
      
      mean.expec= out.optim[[9]]/sum(out.optim[[9]])
      
      h1_optim = sqrt(sum((sqrt(t.dens) - sqrt(optim.pred))^2)) / sqrt(2)
      h2_mean = sqrt(sum((sqrt(t.dens) - sqrt(mean.possum))^2)) / sqrt(2)
      h_3_expec_post = sqrt(sum((sqrt(t.dens) - sqrt(mean.expec))^2)) / sqrt(2)
      
      results <- rbind(results, data.frame(
        method = method[met],
        sample_size = sample.sim[samp],
        iteration = iter,
        h1_optim = h1_optim,
        h2_mean = h2_mean,
        exepec_post = h_3_expec_post
      ))
      cat(iter)
      
    }
  }
}

summary(results$h1_optim)
summary(results$h2_mean)
summary(results$exepec_post)

write.csv(results, file = "/Users/hbolfarine/Dropbox/DSS_MIX/DSS_MIX_06_2024/simulation_results_04_13.csv", row.names = FALSE)

sim_results = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/simulation_results_04_13.csv")

# sim_results = data.frame(results)

library(tidyverse)
library(reshape2)

# sim_tab = sim_results %>% 
#   group_by(method,as.factor(sample_size)) %>% 
#   summarise(mean1 = mean(h1_optim))
# 
# sim_tab = sim_results %>% 
#   group_by(method,as.factor(sample_size)) %>% 
#   summarise(mean2 = mean(as.numeric(h2_mean)))

sel_results = sim_results[,c(2,4,5,6)]

sel_melt = melt(sel_results,id.vars = c("sample_size"))

# Keep both optim, and mean
# sel_results = sim_results %>% 
#   select(-c(iteration))
# 
# sel_melt = melt(sel_results,id.vars = c("method","sample_size"))

# ggplot(sel_melt, aes(x = as.factor(sample_size), y = value, fill = method)) +
#   geom_boxplot() +
#   facet_wrap(~ method + variable) 

# sel_melt = sel_melt

# sel_melt2 = sel_melt2 %>%  
#   mutate(value = (sel_melt$value*sqrt(2))/2)

library(latex2exp)
hellinger = ggplot(sel_melt, aes(x = as.factor(sample_size), y = value, fill = variable)) +
  geom_boxplot() +
  labs(
    x = "Sample Size",                   # Change x-axis label
    y = TeX("Hellinger distance"),                          # Change y-axis label
    fill = "Method"                  # Change legend title
  ) +
  # scale_fill_grey(start = 0.8, end = 0.2) +
  theme_bw() +
  theme(legend.position = 'bottom') +
  # scale_fill_discrete(labels = c("Summary", "Posterior summary average")) +
  scale_fill_brewer(palette="Greys",labels = c("Posterior Summary Estimate", "Posterior Summary Average", "Posterior Expectation")) 

ggsave("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/paper_test_overleaf/hellinger_DPM_04_14.pdf", 
       plot = hellinger, width = 7, height = 4)

  # geom_jitter(shape=16, position=position_jitter(0.2)) +
  # facet_wrap(~method) 


# ggplot(sim_results, aes(x = as.factor(sample_size), y = h2_mean)) +
#   geom_boxplot() +
#   geom_jitter(shape=16, position=position_jitter(0.2)) +
#   facet_wrap(~method) 



