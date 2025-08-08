library(tidyverse)
library(latex2exp)
library(cowplot)
# library(LaplacesDemon)

dlaplace <- function(x, mu, b){
  1/(2*b) * exp(-abs(x-mu)/b)
}

mixt.dss = function(x,list.param,sel.K){
  sum.final = 0
  sum.temp  = 0 
  for(i in 1:sel.K){
    length.sigma = length(list.param$sigmasq[[sel.K]])
    if(length.sigma == 1){
      sigma.i = 1
    }else{
      sigma.i = i
    }
    sum.temp = list.param$pro[[sel.K]][i]*dnorm(x,mean = list.param$mean[[sel.K]][i],sd = sqrt(list.param$sigmasq[[sel.K]][sigma.i]))
    sum.final = sum.final + sum.temp
  }
  return(sum.final)
}

mixtures_true = function(x, index.possum = ""){
  if(index.possum == "example_01"){
    dens1 =  0.5*dnorm(x,mean = 2,sd = 1) + 0.5*dnorm(x,mean = 6,sd = 1)
    return(dens1)
  }else if(index.possum == "example_02"){
    dens2 =  0.2*dnorm(x,mean = 19,sd = sqrt(5)) + 0.2*dnorm(x,mean = 19,sd = 1) +
      0.25*dnorm(x,mean = 23, sd = 1) + 0.2*dnorm(x,mean = 29, sd = 1) + 0.15*dnorm(x,mean = 33, sd = sqrt(2))
    return(dens2)
  }else if(index.possum == "examp_laplace"){
    # dens3 = 0.5*LaplacesDemon::dlaplace(x,location = 0,scale = 1) + 0.5*LaplacesDemon::dlaplace(x,location = 8, scale  = 1)
    dens3 = 0.4*LaplacesDemon::dlaplace(x,location = -5,scale = 1.5) + 0.6*LaplacesDemon::dlaplace(x,location = 5, scale  = 1)
    return(dens3)
  }else if(index.possum == 4){
    0.75*dnorm(x, mean = 0, sd = 1) + 0.25*dnorm(x, mean = 3/2, sd = 1/3)
  }
}

# Plot densities 

# out.mix = possum_OFMM
# out.dpm = possum_DPM
# out.mfm = possum_MFM

plot_possum_results = function(out.mix,out.dpm,out.mfm, y.data, index = 1, sel.K = 1, plot_density = FALSE){
  
  data.bind = rbind(out.mix[[1]],out.dpm[[1]],out.mfm[[1]])
  
  data.bind$method = factor(data.bind$method,levels = c("OBMM","MFM","DPM"))
  
  p1 = ggplot(data.bind, aes(x = num.fact, y = avg.possum, color = method)) + 
    # geom_hline(yintercept = -0.05, color = "gray") +
    # geom_hline(yintercept = 0.05, color = "gray") +
    geom_vline(xintercept = sel.K, colour="gray"  , linetype = "dashed") +
    geom_hline(yintercept = 0, colour="gray") +
    geom_pointrange(aes(ymin = lower, ymax = upper, shape = method, linetype = method), size = 0.4) + 
    geom_line(aes(x = num.fact, y = avg.possum, linetype = method)) +
    scale_x_continuous(breaks = 1:kmax) +
    scale_shape(solid = FALSE) +
    theme_minimal() +
    theme(legend.position="bottom")
  
  p2 = ggplot(data.bind, aes(num.fact, avg.possum, color = method)) + geom_point(aes(shape = method)) + 
    # geom_hline(yintercept =  0.05, color = "gray") +
    # geom_hline(yintercept = -0.05, color = "gray") +
    geom_vline(xintercept = sel.K, colour="gray"  , linetype = "dashed") +
    geom_hline(yintercept = 0, colour="gray", linetype = "dashed") +
    geom_line(aes(x = num.fact, y = avg.possum, linetype = method)) +
    scale_x_continuous(breaks = 1:kmax) +
    theme_minimal() +
    scale_shape(solid = FALSE) +
    theme(legend.position="bottom")
  # scale_colour_grey()
  
  # Plot densities 
  if(plot_density == TRUE){
    
    min.y = min(y.data)
    max.y = max(y.data)
    
    data.distib = data.frame(y.seq =  seq(min.y, max.y, length.out = 500))
    
    data.distib = data.distib %>% 
      mutate(True = mixtures_true(y.seq,index),
             OBMM = mixt.dss(y.seq,list.param = out.mix[[2]],sel.K = sel.K),
             MFM = mixt.dss(y.seq,list.param = out.mfm[[2]],sel.K = sel.K),
             DPM = mixt.dss(y.seq,list.param = out.dpm[[2]],sel.K = sel.K))
    
    data.distib = data.distib %>% 
      as_tibble() %>% 
      pivot_longer(-1)
    
    data.distib$name = factor(data.distib$name,levels = c("True","OBMM","MFM","DPM"))
    
    # y.data = data.frame(y.data)
    p3 = ggplot(data.distib) +
      geom_histogram(data = y.data, aes(x = y.data, y = ..density..), bins = 30 ,colour="darkgrey", fill="white") +
      geom_line(aes(x = y.seq, y = value, color = name, linetype=name)) +
      scale_colour_manual("",values=c(1,2,3,4)) +
      scale_linetype_manual("",values=c(1,2,3,4)) +
      theme_classic() +
      theme(legend.position="bottom")
    
    p =  (p1 + p2) / p3
    return(p)
    
  }else if(plot_density == "applic"){
    
    min.y = min(y.data)
    max.y = max(y.data)
    
    min.y = min.y-0.15*min.y
    max.y = max.y+0.15*max.y
    
    data.distib = data.frame(y.seq =  seq(min.y, max.y, length.out = 500))
    
    data.distib = data.distib %>% 
      mutate(OBMM = mixt.dss(y.seq,list.param = out.mix[[2]],sel.K = sel.K),
             MFM = mixt.dss(y.seq,list.param = out.mfm[[2]],sel.K = sel.K),
             DPM = mixt.dss(y.seq,list.param = out.dpm[[2]],sel.K = sel.K))
    
    data.distib = data.distib %>% 
      as_tibble() %>% 
      pivot_longer(-1)
    
    data.distib$name = factor(data.distib$name,levels = c("OBMM","MFM","DPM"))
    
    y.data = data.frame(y.data)
    p3 = ggplot(data.distib) +
      geom_histogram(data = y.data, aes(x = y.data, y = ..density..), bins = 35,colour="darkgrey", fill="white") +
      geom_line(aes(x = y.seq, y = value, color = name, linetype=name)) +
      scale_colour_manual("",values=c(2,3,4)) +
      scale_linetype_manual("",values=c(2,3,4)) +
      xlim(min.y-0.1*min.y, max.y+0.1*max.y) +
      theme_classic() +
      theme(legend.position="bottom")
    
    p =  (p1 + p2) / p3
    return(p)  
    
  }else {  
    
    p = p1 + p2  
    return(p)  
    
  }
}

#----------------------------------------##----------------------------------------#
# Plots of the posterior number of components - univariate
#----------------------------------------##----------------------------------------#
# possum_BP = BP.galaxy
# possum_DPM = DPM.galaxy
# possum_MFM = MFM.galaxy
plot.postcomp.uni = function(possum_DPM, possum_MFM, possum_BP, K.true = FALSE){
  
  number_of_comp.DPM = data.frame(n.comp = as.numeric(possum_DPM[[3]]),method = "DPM")
  number_of_comp.1 = number_of_comp.DPM %>% 
    group_by(n.comp,method) %>% 
    summarize(count = n()/dim(number_of_comp.DPM)[1])
  
  number_of_comp.MFM = data.frame(n.comp = as.numeric(possum_MFM[[3]]),method = "MFM")
  number_of_comp.2 = number_of_comp.MFM %>% 
    group_by(n.comp,method) %>% 
    summarize(count = n()/dim(number_of_comp.MFM)[1])
  
  number_of_comp.BP = data.frame(n.comp = as.numeric(possum_BP[[3]]),method = "BP")
  number_of_comp.3 = number_of_comp.BP %>% 
    group_by(n.comp,method) %>% 
    summarize(count = n()/dim(number_of_comp.BP)[1])
  
  number_of_comp = rbind(number_of_comp.1,number_of_comp.2,number_of_comp.3)
  
  max_comp = max(number_of_comp$n.comp)
  
  # number_of_comp = number_of_comp %>%
  #   group_by(method) %>%
  #   complete(n.comp = 1:max_comp, fill = list(count = 0))
  
  number_of_comp.1 = number_of_comp %>%
    filter(method %in% c("MFM","DPM")) %>%
    group_by(method) %>%
    complete(n.comp = 1:15, fill = list(count = 0))
  
  number_of_comp.2 = number_of_comp %>%
    filter(method == "BP") %>%
    group_by(method) %>%
    complete(n.comp = 1:max_comp, fill = list(count = 0))
  
  number_of_comp = rbind(number_of_comp.1,number_of_comp.2)
  
  p =  ggplot(number_of_comp, aes(x = n.comp, y = count, shape = method, linetype = method)) +
    # p =  ggplot(number_of_comp, aes(x = n.comp, y = count, shape = method, linetype = method, col = method)) +
    # p =  ggplot(number_of_comp, aes(x = n.comp, y = count, shape = method, color = method)) +
    # p =  ggplot(number_of_comp, aes(x = n.comp, y = count, color = method)) +
    # geom_point(aes(x = n.comp, y = count, shape = method, col = method)) +
    # geom_line(aes(x = n.comp, y = count, linetype = method, col = method)) +
    geom_point(size = 2) +
    geom_line() +
    # geom_bar(stat = "identity") +
    # scale_x_continuous(breaks = 1:max_comp) +
    # theme(legend.position="bottom") +
    # theme_light() +
    theme_classic() +
    xlab("Number of components") +
    ylab("Percentage (Counts)") +
    ggtitle("Posterior on the number of components") +
    theme(legend.position = c(0.8, 0.8), plot.title = element_text(hjust = 0.5))
  # theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) 
  # scale_colour_grey() +
  # theme(legend.position="bottom") 
  # facet_wrap(~method, scales = "free_x")
  
  if(K.true != FALSE){
    p = p + geom_vline(xintercept = K.true, colour="gray"  , linetype = "dashed")
  }
  
  return(p)
  # x = "Number of Components",
  # y = "Percentage (Counts)"
}

#----------------------------------------##----------------------------------------#
# Plots of the posterior number of components - multivariate
#----------------------------------------##----------------------------------------#
# possum_DPM_mult = temp.DPM.mult
# possum_MFM_mult = temp.MFM.mult
# possum_SFM_mult = temp.SFM.mult

plot.postcomp.mult = function(possum_DPM_mult, possum_MFM_mult, possum_SFM_mult, K.true = FALSE, kmax = 10){
  
  number_of_comp.DPM = data.frame(n.comp = as.numeric(unlist(possum_DPM_mult[[3]])),method = "DPM")
  number_of_comp.1 = number_of_comp.DPM %>% 
    group_by(n.comp,method) %>% 
    summarize(count = n()/dim(number_of_comp.DPM)[1])
  
  number_of_comp.MFM = data.frame(n.comp = as.numeric(unlist(possum_MFM_mult[[3]])),method = "MFM")
  number_of_comp.2 = number_of_comp.MFM %>% 
    group_by(n.comp,method) %>% 
    summarize(count = n()/dim(number_of_comp.MFM)[1])
  
  number_of_comp.SFM = data.frame(n.comp = as.numeric(unlist(possum_SFM_mult[[3]])),method = "SFM")
  number_of_comp.3 = number_of_comp.SFM %>% 
    group_by(n.comp,method) %>% 
    summarize(count = n()/dim(number_of_comp.SFM)[1])
  
  number_of_comp = rbind(number_of_comp.1,number_of_comp.2,number_of_comp.3)
  
  max_comp_DPM = max(number_of_comp[which(number_of_comp$method =="DPM"),]$n.comp)
  max_comp_MFM = max(number_of_comp[which(number_of_comp$method =="MFM"),]$n.comp)
  max_comp_SFM = max(number_of_comp[which(number_of_comp$method =="SFM"),]$n.comp)
  
  # number_of_comp = number_of_comp %>%
  #   group_by(method) %>%
  #   complete(n.comp = 1:max_comp, fill = list(count = 0))
  
  number_of_comp_DPM = number_of_comp %>%
    group_by(method) %>% 
    filter(method == "DPM") %>% 
    complete(n.comp = 1:10, fill = list(count = 0))
  
  max_comp_MFM = number_of_comp %>%
    group_by(method) %>% 
    filter(method == "MFM") %>% 
    complete(n.comp = 1:10, fill = list(count = 0))
  
  max_comp_SFM = number_of_comp %>%
    group_by(method) %>% 
    filter(method == "SFM") %>% 
    complete(n.comp = 1:10, fill = list(count = 0))
  
  number_of_comp_all = rbind(number_of_comp_DPM,max_comp_MFM,max_comp_SFM)
  
  # x = "Number of Components",
  # y = "Percentage (Counts)"
  
  p =  ggplot(number_of_comp_all) +
    # geom_point(aes(x = n.comp, y = count, shape = method, col = method)) +
    # geom_line(aes(x = n.comp, y = count, linetype = method, col = method)) +
    geom_point(aes(x = n.comp, y = count, shape = method), size = 2) +
    geom_line(aes(x = n.comp, y = count, linetype = method)) +
    scale_x_continuous(breaks = 1:max(number_of_comp_all$n.comp)) +
    theme(legend.position="bottom") +
    theme_classic() +
    # theme_light() +
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
    # theme(legend.position = c(0.85, 0.85)) #+
    ggtitle("Posterior on the number of components") +
    xlab("Number of components") +
    ylab("Percentage (Counts)") +
    theme(legend.position = c(0.8, 0.8), plot.title = element_text(size = 12, hjust = 0.5), 
          axis.text.x = element_text(size = 8))
  # theme(legend.position="bottom") +
  # facet_wrap(~method, scales = "free_x")
  
  if(K.true != FALSE){
    p = p + geom_vline(xintercept = K.true, colour="gray"  , linetype = "dashed")
  }
  
  return(p)
  
}

#----------------------------------------##----------------------------------------#
# Plot multivariate results
#----------------------------------------##----------------------------------------#

# out.sfm.mult = temp.SFM.mult
# out.dpm.mult = temp.DPM.mult
# out.mfm.mult = temp.MFM.mult

plot.estim.comp.mult = function(out.sfm.mult,out.dpm.mult,out.mfm.mult, y.data, index = 1, sel.K = 1, plot_density = FALSE, title = "", kmax = 10){
  
  data.bind = rbind(out.dpm.mult[[1]],out.mfm.mult[[1]],out.sfm.mult[[1]])
  data.bind$method = factor(data.bind$method,levels = c("DPM.mult","MFM.mult","SFM.mult"))
  
  data.bind  = data.bind %>% 
    mutate(method = fct_recode(method, "DPM" = "DPM.mult", "MFM" = "MFM.mult", "SFM" = "SFM.mult"))
  
  # p = ggplot(data.bind, aes(x = num.fact, y = avg.possum, color = method)) + 
  #   # geom_hline(yintercept = -0.05, color = "gray") +
  #   # geom_hline(yintercept = 0.05, color = "gray") +
  #   geom_vline(xintercept = sel.K, colour="gray", linetype = "dashed") +
  #   geom_hline(yintercept = 0, colour="gray") +
  #   scale_x_continuous(breaks = 1:kmax) +
  #   geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.3, position = position_dodge(0.6)) +
  #   theme_bw() +
  #   theme(legend.position="bottom") +
  #   theme(legend.position="bottom", legend.title=element_blank()) +
  #   ggtitle(title) +
  #   xlab("Number of components") +
  #   ylab("func in latex") +
  #   # geom_line(aes(x = num.fact, y = avg.possum, linetype = method)) +
  #   # scale_colour_manual(values=c("red","green","blue")) +
  #   scale_x_continuous(breaks = 1:kmax) +
  #   scale_linetype_manual(values=c(2,3,4)) +
  #   scale_shape(solid = FALSE) +
  #   theme(legend.position="bottom") +
  #   theme(legend.position="bottom", legend.title=element_blank()) +
  #   theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
  #   facet_wrap(~ method)
  
  # p = ggplot(data.bind, aes(x = num.fact, y = avg.possum, color = method)) +
  p = ggplot(data.bind, aes(x = num.fact, y = avg.possum)) +
    # p = ggplot(data.bind, aes(x = num.fact, y = avg.possum)) + 
    # geom_hline(yintercept = -0.05, color = "gray") +
    geom_hline(yintercept = 0, colour="gray", linetype = "dashed") +
    scale_x_continuous(breaks = 1:kmax) +
    # geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.2, position = position_dodge(0.6)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.4) +
    geom_point(aes(x = num.fact, y = avg.possum),size = 2.5, shape = 1) +
    # theme(legend.position="bottom") +) +
    theme(legend.position="bottom") +
    theme(legend.position="bottom", legend.title=element_blank()) +
    # ggtitle(title) +
    ggtitle("Posterior summary discrepancy function") +
    xlab(TeX("Surrogate dimension $k$")) +
    ylab(TeX("$d_{n}^{k}$")) +
    # geom_line(aes(x = num.fact, y = avg.possum, linetype = method)) +
    geom_line(aes(x = num.fact, y = avg.possum)) +
    # scale_colour_manual(values=c("red","green","blue")) +
    scale_x_continuous(breaks = 1:kmax) +
    scale_linetype_manual(values=c(2,3,4)) +
    scale_shape(solid = FALSE) +
    # theme_light() +
    # theme(legend.position="bottom", legend.title=element_blank()) +
    # theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.position="bottom", legend.title=element_blank(),
          axis.text.x = element_text(size = 8),  # Rotate x-axis labels
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 12, hjust = 0.5), 
          plot.margin = unit(c(.2,.2,.2,.5), "cm")) + # margin(t, r, b, l )
    # theme_classic() +
    # theme_light() +
    # theme_minimal() +
    # theme(legend.position="bottom") +
    facet_wrap(~method)
  
  # p2 = ggplot(data.bind, aes(num.fact, avg.possum, color = method)) + geom_point(aes(shape = method)) +
  #   geom_hline(yintercept =  0.05, color = "gray") +
  #   geom_hline(yintercept = -0.05, color = "gray") +
  #   geom_vline(xintercept = sel.K, colour="gray"  , linetype = "dashed") +
  #   geom_hline(yintercept = 0, colour="gray", linetype = "dashed") +
  #   geom_line(aes(x = num.fact, y = avg.possum, linetype = method)) +
  #   scale_x_continuous(breaks = 1:kmax) +
  #   theme_classic() +
  #   scale_shape(solid = FALSE) +
  #   theme(legend.position="bottom")
  # # scale_colour_grey()
  
  # # Plot densities 
  # if(plot_density == "applic"){
  #   cat("applic")
  #   if(index == 3){  
  #     # Faithful
  #     x.dist <- seq(1, 7, length.out = 100)
  #     y.dist <- seq(35, 100, length.out = 100)
  #     dat <- expand_grid(x = x.dist, y = y.dist)
  #   }
  #   #OBMM
  #   z2 = numeric(dim(dat)[1])
  #   z3 = numeric(dim(dat)[1])
  #   z4 = numeric(dim(dat)[1])
  #   for(i in 1:dim(dat)[1]){
  #     z2[i] = mixt.dss.mult(x = dat[i,1], y = dat[i,2],list.param = out.mix[[2]],sel.K)
  #     z3[i] = mixt.dss.mult(x = dat[i,1], y = dat[i,2],list.param = out.mfm.mult[[2]],sel.K)
  #     z4[i] = mixt.dss.mult(x = dat[i,1], y = dat[i,2],list.param = out.dpm.mult[[2]],sel.K)
  #   }
  #   
  #   data_distrib = data.frame(dat, OFMM_mult = z2, MFM_mult = z3, DPM.mult = z4)
  #   
  #   data_distrib = data_distrib %>% 
  #     as_tibble() %>% 
  #     pivot_longer(-c(1,2))
  #   
  #   y.data = data.frame(X = y.data[,1], Y = y.data[,2])
  #   
  #   data_distrib$name = factor(data_distrib$name,levels = c("OFMM_mult","MFM_mult","DPM.mult"))
  #   
  #   p3 = ggplot(data_distrib) + 
  #     geom_point(data = y.data, aes(x = X,y = Y),col = "darkgray",shape = 3, size = 1) +
  #     stat_contour(aes(x = x, y = y, z = value, color = name, linetype=name)) +
  #     theme_classic()
  #   
  #   p =  (p1 + p2) / p3
  #   return(p)
  #   
  # }else if(plot_density == "simulation"){
  #   cat("simulation")
  #   
  #   if(index == 1){
  #     # Simulation 1
  #     x.dist <- seq(-3, 15, length.out = 200)
  #     y.dist <- seq(-3, 20, length.out = 200)
  #     dat <- expand_grid(x = x.dist, y = y.dist)
  #     
  #     z1 = numeric(dim(dat)[1])
  #     for(i in 1:dim(dat)[1]){
  #       z1[i] = mixtures_true_mult(x = dat[i,1], y = dat[i,2], index = index)
  #     }
  #     
  #   }else if(index == 2){
  #     # Simulation 2
  #     x.dist <- seq(-2, 10, length.out = 200)
  #     y.dist <- seq(-2, 20, length.out = 200)
  #     dat <- expand_grid(x = x.dist, y = y.dist)
  #     
  #     z1 = numeric(dim(dat)[1])
  #     for(i in 1:dim(dat)[1]){
  #       z1[i] = mixtures_true_mult(x = dat[i,1], y = dat[i,2], index = index)
  #     }
  #   }
  #   
  #   z2 = numeric(dim(dat)[1])
  #   z3 = numeric(dim(dat)[1])
  #   z4 = numeric(dim(dat)[1])
  #   for(i in 1:dim(dat)[1]){
  #     z2[i] = mixt.dss.mult(x = dat[i,1], y = dat[i,2],list.param = out.mix[[2]],sel.K)
  #     z3[i] = mixt.dss.mult(x = dat[i,1], y = dat[i,2],list.param = out.mfm.mult[[2]],sel.K)
  #     z4[i] = mixt.dss.mult(x = dat[i,1], y = dat[i,2],list.param = out.dpm.mult[[2]],sel.K)
  #   }  
  #   
  #   data_distrib = data.frame(dat,True = z1, OFMM_mult = z2, MFM_mult = z3, DPM.mult = z4)
  #   
  #   data_distrib = data_distrib %>% 
  #     as_tibble() %>% 
  #     pivot_longer(-c(1,2))
  #   
  #   data_distrib$name = factor(data_distrib$name,levels = c("True","OFMM_mult","MFM_mult","DPM.mult"))
  #   
  #   y.data = data.frame(X = y.data[,1], Y = y.data[,2])
  #   
  #   p3 = ggplot(data_distrib) + 
  #     geom_point(data = y.data, aes(x = X,y = Y),col = "darkgray",shape = 3, size = 1) +
  #     stat_contour(aes(x = x, y = y, z = value, color = name, linetype=name)) +
  #     theme_classic()
  #   
  #   cat("here")
  #   p =  (p1 + p2) / p3
  #   return(p)
  # }
  # 
  # if(plot_density =="FALSE"){
  #   cat("here")
  #   p =  p1 + p2
  #   return(p)
  # }
  if(sel.K != 0){
  p = p + geom_vline(xintercept = sel.K, colour="gray", linetype = "dashed") 
  }else if(sel.K == 0){
    p = p
  }  
  return(p)
}

#----------------------------------------##----------------------------------------#
# Plot discrepancy function for univariate model
#----------------------------------------##----------------------------------------#

plot.possum.comp.uni = function(out.BP,out.DPM,out.MFM, sel.K = 1, kmax = 10){
  
  data.bind = rbind(out.BP[[1]],out.DPM[[1]],out.MFM[[1]])
  data.bind$method = factor(data.bind$method,levels = c("BP","DPM","MFM"))
  
  # data.bind  = data.bind %>% 
  #   mutate(method = fct_recode(method, "DPM" = "DPM.mult", "MFM" = "MFM.mult", "SFM" = "SFM.mult"))
  # y.data = data.frame(y.data)
  # p = ggplot(data.bind, aes(x = num.fact, y = avg.possum, color = method)) +
  p = ggplot(data.bind, aes(x = num.fact, y = avg.possum)) +
    # geom_hline(yintercept = -0.05, color = "gray") +
    # geom_hline(yintercept = 0.05, color = "gray") +
    geom_vline(xintercept = sel.K, colour="gray", linetype = "dashed") +
    geom_hline(yintercept = 0, colour="gray", linetype = "dashed") +
    scale_x_continuous(breaks = 1:kmax) +
    # geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.2, position = position_dodge(0.6)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5) +
    geom_point(aes(x = num.fact, y = avg.possum),size = 2.5, shape = 1) +
    theme(legend.position="bottom") +
    theme(legend.position="bottom", legend.title=element_blank()) +
    ggtitle("Posterior summary discrepancy function") +
    xlab(TeX("Surrogate dimension $k$")) +
    ylab(TeX("$d_{n}^{k}$")) +
    # ylab(TeX("Formula: $\\frac{2hc^2}{\\delta^\\beta}$")) +
    # geom_line(aes(x = num.fact, y = avg.possum, linetype = method)) +
    geom_line(aes(x = num.fact, y = avg.possum)) +
    # scale_colour_manual(values=c("red","green","blue")) +
    scale_x_continuous(breaks = 1:kmax) +
    scale_linetype_manual(values=c(2,3,4)) +
    scale_shape(solid = FALSE) +
    # theme_light() +
    # theme(legend.position="bottom", legend.title=element_blank()) +
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.position="bottom", legend.title=element_blank(),
          axis.text.x = element_text(size = 10),  # Rotate x-axis labels
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 12, hjust = 0.5),
          plot.margin = unit(c(.2,.2,.2,.5), "cm")) + # margin(t, r, b, l ) +
    # theme_classic() +
    # theme_minimal() +
    # theme(legend.position="bottom") +
    facet_wrap(~method)
  
  return(p)
}

#----------------------------------------##----------------------------------------#
# Plot discrepancy function for univariate model
#----------------------------------------##----------------------------------------#

plot.possum.uni = function(out, sel.K = FALSE, kmax = 10, method = "", y.lim = c(-2,0.5)){
  
  # data.bind  = data.bind %>% 
  #   mutate(method = fct_recode(method, "DPM" = "DPM.mult", "MFM" = "MFM.mult", "SFM" = "SFM.mult"))
  # y.data = data.frame(y.data)
  # p = ggplot(data.bind, aes(x = num.fact, y = avg.possum, color = method)) +
  p = ggplot(out, aes(x = num.fact, y = avg.possum)) +
    # geom_hline(yintercept = -0.05, color = "gray") +
    # geom_hline(yintercept = 0.05, color = "gray") +
    # geom_vline(xintercept = sel.K, colour="gray", linetype = "dashed") +
    geom_hline(yintercept = 0, colour="lightgray", linetype = "dashed") +
    scale_x_continuous(breaks = 1:kmax) +
    # geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.2, position = position_dodge(0.6)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.4) +
    geom_point(aes(x = num.fact, y = avg.possum),size = 2.5, shape = 1) +
    # geom_point(aes(x = num.fact, y = avg.possum),size = 2.0) +
    # theme(legend.position="bottom") +
    ggtitle("Posterior summary discrepancy function") +
    xlab(TeX("Summary dimension $k$")) +
    ylab(TeX("$d_{n}^{k}$")) +
    # ylab(TeX("Formula: $\\frac{2hc^2}{\\delta^\\beta}$")) +
    # geom_line(aes(x = num.fact, y = avg.possum, linetype = method)) +
    geom_line(aes(x = num.fact, y = avg.possum)) +
    # scale_colour_manual(values=c("red","green","blue")) +
    scale_x_continuous(breaks = 1:kmax) +
    scale_linetype_manual(values=c(2,3,4)) +
    scale_shape(solid = FALSE) +
    # theme_light() +
    # theme(legend.position="bottom", legend.title=element_blank()) +
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
    # theme_bw() #+
    theme_bw() +
    ylim(y.lim[1], y.lim[2]) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.position="bottom", legend.title=element_blank(),
          axis.text.x = element_text(size = 10),  # Rotate x-axis labels
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 11, hjust = 0.5),
          plot.margin = unit(c(.2,.2,.2,.5), "cm") # margin(t, r, b, l )
    )
  # theme_minimal() #+
  # theme(legend.position="bottom") +
  # facet_wrap(~method)
  
  if(sel.K != FALSE){
    p = p + geom_vline(xintercept = sel.K, colour="gray", linetype = "dashed")
  }
  
  return(p)
}

#----------------------------------------##----------------------------------------#
# Plot point estimate prediction densities for univariate model
#----------------------------------------##----------------------------------------#

# out.MFM = temp.MFM
# out.BP = temp.BP
# out.DPM = temp.DPM

plot.estim.uni = function(out.BP,out.DPM,out.MFM, sel.K = 1, title = "", kmax = 10, index = "", plot_density = TRUE, y.data){
  
  data.bind = rbind(out.BP[[1]],out.DPM[[1]],out.MFM[[1]])
  data.bind$method = factor(data.bind$method,levels = c("BP","DPM","MFM"))
  
  if(plot_density == TRUE){
    
    min.y = min(y.data)
    max.y = max(y.data)
    
    data.distib = data.frame(y.seq =  seq(min.y, max.y, length.out = 500))
    
    data.distib = data.distib %>% 
      mutate(True = mixtures_true(y.seq,index),
             BP = mixt.dss(y.seq,list.param = out.BP[[2]],sel.K = sel.K),
             MFM = mixt.dss(y.seq,list.param = out.MFM[[2]],sel.K = sel.K),
             DPM = mixt.dss(y.seq,list.param = out.DPM[[2]],sel.K = sel.K))
    
    data.distib = data.distib %>% 
      as_tibble() %>% 
      pivot_longer(-1)
    
    data.distib$name = factor(data.distib$name,levels = c("True","BP","MFM","DPM"))
    
    y.data = data.frame(y.data)
    p1 = ggplot(data.distib) +
      geom_histogram(data = y.data, aes(x = y.data, y = ..density..), bins = 30 ,colour="darkgrey", fill="white") +
      geom_line(aes(x = y.seq, y = value, color = name, linetype=name)) +
      scale_colour_manual("",values=c(1,2,3,4)) +
      scale_linetype_manual("",values=c(1,2,3,4)) +
      theme_classic() +
      theme(legend.position="bottom")
    
    p =  p1
    
  }else if(plot_density == "applic"){
    
    min.y = min(y.data)
    max.y = max(y.data)
    
    min.y = min.y-0.15*min.y
    max.y = max.y+0.15*max.y
    
    data.distib = data.frame(y.seq =  seq(min.y, max.y, length.out = 500))
    
    data.distib = data.distib %>% 
      mutate(BP = mixt.dss(y.seq,list.param = out.BP[[2]],sel.K = sel.K),
             MFM = mixt.dss(y.seq,list.param = out.MFM[[2]],sel.K = sel.K),
             DPM = mixt.dss(y.seq,list.param = out.DPM[[2]],sel.K = sel.K))
    
    data.distib = data.distib %>% 
      as_tibble() %>% 
      pivot_longer(-1)
    
    data.distib$name = factor(data.distib$name,levels = c("BP","MFM","DPM"))
    
    y.data = data.frame(y.data)
    p2 = ggplot(data.distib) +
      geom_histogram(data = y.data, aes(x = y.data, y = ..density..), bins = 35,colour="darkgrey", fill="white") +
      geom_line(aes(x = y.seq, y = value, color = name, linetype=name)) +
      scale_colour_manual("",values=c(2,3,4)) +
      scale_linetype_manual("",values=c(2,3,4)) +
      xlim(min.y-0.1*min.y, max.y+0.1*max.y) +
      theme_classic() +
      theme(legend.position="bottom")
    
    p = p2
    return(p)  
    
  }
  
}

#----------------------------------------##----------------------------------------#
# Plot point estimate prediction densities for univariate model
#----------------------------------------##----------------------------------------#

dc.possum.uncertain.quant = function(possum_BP, possum_MFM, possum_DPM, K.sel, y.data, 
                                     index = FALSE, param.summ = TRUE, split.method = FALSE,
                                     scale.plot = FALSE, quant = c(0.025,0.975)){
  
  y.fit.BP = possum_BP[[4]]
  y.fit.MFM  = possum_MFM[[4]]
  y.fit.DPM = possum_DPM[[4]]
  
  N = ncol(y.fit.BP)
  H = nrow(y.fit.BP)
  
  pred.fit.BP = pred.fit.MFM = pred.fit.DPM = matrix(0,N,H)
  mean_post.BP = sigmasq_post.BP= prob_post.BP = matrix(0,N,K.sel)
  mean_post.MFM = sigmasq_post.MFM = prob_post.MFM= matrix(0,N,K.sel)
  mean_post.DPM = sigmasq_post.DPM = prob_post.DPM = matrix(0,N,K.sel)
  
  # BP
  for(i in 1:N){
    dens <- densityMclust(y.fit.BP[i,], G = K.sel,plot = F, verbose = F)
    dens1 = dens$density
    
    # Save parameters for plot.
    prob_post.BP[i,] = dens$parameters$pro
    mean_post.BP[i,] = dens$parameters$mean
    sigmasq_post.BP[i,] = dens$parameters$variance$sigmasq
    pred.fit.BP[i,] = dens1
  }
  
  # MFM
  for(i in 1:N){
    dens <- densityMclust(y.fit.MFM[i,], G = K.sel,plot = F, verbose = F)
    dens1 = dens$density
    
    # Save parameters for plot.
    prob_post.MFM[i,] = dens$parameters$pro
    mean_post.MFM[i,] = dens$parameters$mean
    sigmasq_post.MFM[i,] = dens$parameters$variance$sigmasq
    pred.fit.MFM[i,] = dens1
  }
  
  # DPM
  for(i in 1:N){
    dens <- densityMclust(y.fit.DPM[i,], G = K.sel,plot = F, verbose = F)
    dens1 = dens$density
    
    # Save parameters for plot.
    prob_post.DPM[i,] = dens$parameters$pro
    mean_post.DPM[i,] = dens$parameters$mean
    sigmasq_post.DPM[i,] = dens$parameters$variance$sigmasq
    pred.fit.DPM[i,] = dens1
  }
  
  mu.seq = paste0("mu", 1:K.sel)
  sigma.seq = paste0("sigma", 1:K.sel)
  prob.seq = paste0("prob", 1:K.sel)
  
  #------------------------------------#
  
  mean_post.BP = data.frame(mean_post.BP)
  colnames(mean_post.BP) = mu.seq
  mean_est_post.BP = apply(mean_post.BP,2, mean)
  
  mean_post.MFM = data.frame(mean_post.MFM)
  colnames(mean_post.MFM) = mu.seq
  mean_est_post.MFM = apply(mean_post.MFM,2, mean)
  
  mean_post.DPM = data.frame(mean_post.DPM)
  colnames(mean_post.DPM) = mu.seq
  mean_est_post.DPM = apply(mean_post.DPM,2, mean)
  
  #------------------------------------#
  
  sigmasq_post.BP = data.frame(sigmasq_post.BP)
  colnames(sigmasq_post.BP) = sigma.seq
  sigmasq_est_post.BP = apply(sigmasq_post.BP,2, mean)
  
  sigmasq_post.MFM = data.frame(sigmasq_post.MFM)
  colnames(sigmasq_post.MFM) = sigma.seq
  sigmasq_est_post.MFM = apply(sigmasq_post.MFM,2, mean)
  
  sigmasq_post.DPM = data.frame(sigmasq_post.DPM)
  colnames(sigmasq_post.DPM) = sigma.seq
  sigmasq_est_post.DPM = apply(sigmasq_post.DPM,2, mean)
  
  #------------------------------------#
  
  prob_post.BP = data.frame(prob_post.BP)
  colnames(prob_post.BP) = prob.seq
  prob_est_post.BP = apply(prob_post.BP,2, mean)
  
  prob_post.MFM = data.frame(prob_post.MFM)
  colnames(prob_post.MFM) = prob.seq
  prob_est_post.MFM = apply(prob_post.MFM,2, mean)
  
  prob_post.DPM = data.frame(prob_post.DPM)
  colnames(prob_post.DPM) = prob.seq
  prob_est_post.DPM = apply(prob_post.DPM,2, mean)
  
  #------------------------------------#
  #------------------------------------#
  
  mean.post.BP = mean_post.BP %>% 
    mutate(method = "BP")
  
  mean.post.MFM = mean_post.MFM %>% 
    mutate(method = "MFM")
  
  mean.post.DPM = mean_post.DPM %>% 
    mutate(method = "DPM")
  
  #------------------------------------#
  
  sigmasq.post.BP = sigmasq_post.BP %>% 
    mutate(method = "BP")
  
  sigmasq.post.MFM = sigmasq_post.MFM %>% 
    mutate(method = "MFM")
  
  sigmasq.post.DPM = sigmasq_post.DPM %>% 
    mutate(method = "DPM")
  
  #------------------------------------#
  
  prob.post.BP = prob_post.BP %>% 
    mutate(method = "BP")
  
  prob.post.MFM = prob_post.MFM %>% 
    mutate(method = "MFM")
  
  prob.post.DPM = prob_post.DPM %>% 
    mutate(method = "DPM")
  
  #------------------------------------#
  #------------------------------------#
  
  mean.melt.OGM = melt(mean.post.BP)
  mean.melt.MFM = melt(mean.post.MFM)
  mean.melt.DPM = melt(mean.post.DPM)
  
  #------------------------------------#
  
  sigmasq.melt.BP = melt(sigmasq.post.BP)
  sigmasq.melt.MFM = melt(sigmasq.post.MFM)
  sigmasq.melt.DPM = melt(sigmasq.post.DPM)
  
  #------------------------------------#
  
  prob.melt.BP = melt(prob.post.BP)
  prob.melt.MFM = melt(prob.post.MFM)
  prob.melt.DPM = melt(prob.post.DPM)
  
  #------------------------------------#
  #------------------------------------#  
  
  mean.mealt.all = rbind(mean.melt.OGM,mean.melt.MFM,mean.melt.DPM)
  mean.mealt.all = mean.mealt.all %>% 
    mutate(param = "mu")
  var.melt.all = rbind(sigmasq.melt.BP,sigmasq.melt.MFM,sigmasq.melt.DPM)
  var.melt.all = var.melt.all %>% 
    mutate(param = "sigma")
  prob.melt.all = rbind(prob.melt.BP,prob.melt.MFM,prob.melt.DPM)
  prob.melt.all = prob.melt.all %>% 
    mutate(param = "prob")
  
  # dens2 =  0.2*dnorm(x,mean = 16,sd = sqrt(5)) + 0.2*dnorm(x,mean = 19,sd = 1) +
  #   0.25*dnorm(x,mean = 23, sd = 1) + 0.2*dnorm(x,mean = 29, sd = 1) + 0.15*dnorm(x,mean = 33, sd = sqrt(1))
  # 
  
  param.all = rbind(mean.mealt.all,var.melt.all,prob.melt.all)
  
  if(index != FALSE){
    
    if(index == "example2"){
      data.lines.mu = data.frame(vec = c(19,19,23,29,33), variable = mu.seq, param = "mu")
      data.lines.sigma = data.frame(vec = c(5,1,1,1,2), variable = sigma.seq, param = "sigma")
      data.lines.prob = data.frame(vec = c(0.2,0.2,0.25,0.2,0.15), variable = prob.seq, param = "prob")
    }else if(index == "laplace"){
      data.lines.mu = data.frame(vec = c(0,8), variable = mu.seq, param = "mu")
      data.lines.sigma = data.frame(vec = c(1,1), variable = sigma.seq, param = "sigma")
      data.lines.prob = data.frame(vec = c(0.5,0.5), variable = prob.seq, param = "prob")
    }else if(index == "wand"){
      data.lines.mu = data.frame(vec = c(0,2/3), variable = mu.seq, param = "mu")
      data.lines.sigma = data.frame(vec = c(1,1/3), variable = sigma.seq, param = "sigma")
      data.lines.prob = data.frame(vec = c(0.75,0.25), variable = prob.seq, param = "prob")
    }
    data.lines = rbind(data.lines.mu,data.lines.sigma,data.lines.prob)
  }
  #BP
  data.estm.mu.BP = data.frame(vec = mean_est_post.BP, variable = mu.seq, param = "mu", method = "BP")
  data.estm.sigma.BP = data.frame(vec = sigmasq_est_post.BP, variable = sigma.seq, param = "sigma", method = "BP")
  data.estm.prob.BP = data.frame(vec = prob_est_post.BP, variable = prob.seq, param = "prob", method = "BP")
  
  data.estm.mu.MFM = data.frame(vec = mean_est_post.MFM, variable = mu.seq, param = "mu", method = "MFM")
  data.estm.sigma.MFM = data.frame(vec = sigmasq_est_post.MFM, variable = sigma.seq, param = "sigma", method = "MFM")
  data.estm.prob.MFM = data.frame(vec = prob_est_post.MFM, variable = prob.seq, param = "prob", method = "MFM")
  
  data.estm.mu.DPM = data.frame(vec = mean_est_post.DPM, variable = mu.seq, param = "mu", method = "DPM")
  data.estm.sigma.DPM = data.frame(vec = sigmasq_est_post.DPM, variable = sigma.seq, param = "sigma", method = "DPM")
  data.estm.prob.DPM = data.frame(vec = prob_est_post.DPM, variable = prob.seq, param = "prob", method = "DPM")
  
  data.points = rbind(data.estm.mu.BP,data.estm.sigma.BP,data.estm.prob.BP,
                      data.estm.mu.MFM,data.estm.sigma.MFM,data.estm.prob.MFM,
                      data.estm.mu.DPM,data.estm.sigma.DPM,data.estm.prob.DPM)
  
  if(param.summ == TRUE){
    # p1 = ggplot(param.all) +
    #   geom_density(aes(x = value, fill = method, col = method),alpha = 0.4) +
    #   facet_wrap(param~variable, scales = "free", ncol = K.sel) +
    #   # geom_vline(data = data.lines, aes(xintercept = vec), col = "red") +
    #   labs(title = "Density Plots", x = "Value", y = "Density") +
    #   theme_minimal() +
    #   theme(legend.position="bottom")
    
    # Violin plot 
    # http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization
    
    data_summary <- function(x){
      m <- mean(x)
      # ymin <- m-sd(x)
      # ymax <- m+sd(x)
      ymin <- as.numeric(quantile(x,probs = 0.975))
      ymax <- as.numeric(quantile(x,probs = 0.025))
      return(c(y=m,ymin=ymin,ymax=ymax))
    }
    
    # The bar in the violin plot indicates the 95% credible interval, and
    # the point indicates the posterior mean. 
    
    p2 = ggplot(param.all) +
      # geom_violin(aes(x = method, y = value),alpha = 0.5, draw_quantiles = c(0.025,0.975)) +
      geom_violin(aes(x = method, y = value),alpha = 0.5) +
      facet_wrap(param~variable, scales = "free", ncol = K.sel) +
      labs(title = "Density Plots", x = "Value", y = "Density") +
      theme_minimal() +
      # geom_boxplot(aes(x = method, y = value),width=0.1, outlier.shape = NA) +
      # stat_summary(aes(x = method, y = value), fun.y=mean, geom="point", size=2, color="red", shape = 4) +
      # stat_summary(aes(x = method, y = value),fun.data=mean_sdl, mult=1,  geom="pointrange", color="red")
      stat_summary(aes(x = method, y = value),fun.data=data_summary, color="black", size=0.3, shape = 20) +
      geom_hline(data = data.lines, aes(yintercept = vec), col = "red") +
      # coord_flip() +
      theme(legend.position="bottom")
    
    # p1 = p1 + geom_vline(data = data.points, aes(xintercept = vec, col = method))
    
    if(index != FALSE){
      # p1 = p1 + geom_vline(data = data.lines, aes(xintercept = vec), col = "red")
      p2 = p2 + geom_hline(data = data.lines, aes(yintercept = vec), col = "red")
    }
    
    # p = p1/p2
    p = p2
  }
  
  #------------------------------------#
  #------------------------------------#  
  
  if(param.summ == FALSE){
    
    if(sum(scale.plot) == FALSE){
      min.y = min(y.data) - 0.1*min(y.data)
      max.y = max(y.data) + 0.1*max(y.data)
    }else{
      min.y = scale.plot[1]
      max.y = scale.plot[2]
    }
    y.seq = seq(min.y, max.y, length.out = 500)
    
    if(index != FALSE){
      true.dens = mixtures_true(y.seq,index)
    }  
    
    list.param.BP = list(mean = mean_post.BP, sigmasq = sigmasq_post.BP, pro = prob_post.BP)
    list.param.MFM = list(mean = mean_post.MFM, sigmasq = sigmasq_post.MFM, pro = prob_post.MFM)
    list.param.DPM = list(mean = mean_post.DPM, sigmasq = sigmasq_post.DPM, pro = prob_post.DPM)
    
    N = nrow(list.param.BP$mean)
    predi.mat.BP = matrix(0,N,length(y.seq))
    for(i in 1:N){
      for(j in 1:length(y.seq)){
        predi.mat.BP[i,j] = sum(list.param.BP$pro[i,]*dnorm(y.seq[j],mean = as.numeric(list.param.BP$mean[i,]),sd = as.numeric(sqrt(list.param.BP$sigmasq[i,]))))
      }
      # cat(i)
    }
    
    predi.mat.MFM = matrix(0,N,length(y.seq))
    for(i in 1:N){
      for(j in 1:length(y.seq)){
        predi.mat.MFM[i,j] = sum(list.param.MFM$pro[i,]*dnorm(y.seq[j],mean = as.numeric(list.param.MFM$mean[i,]),sd = as.numeric(sqrt(list.param.MFM$sigmasq[i,]))))
      }
      # cat(i)
    }
    
    predi.mat.DPM = matrix(0,N,length(y.seq))
    for(i in 1:N){
      for(j in 1:length(y.seq)){
        predi.mat.DPM[i,j] = sum(list.param.DPM$pro[i,]*dnorm(y.seq[j],mean = as.numeric(list.param.DPM$mean[i,]),sd = as.numeric(sqrt(list.param.DPM$sigmasq[i,]))))
      }
      # cat(i)
    }
    
    mean.pred.BP = apply(predi.mat.BP, 2, mean)
    mean.pred.MFM = apply(predi.mat.MFM, 2, mean)
    mean.pred.DPM = apply(predi.mat.DPM, 2, mean)
    
    # Apply the function to each column of the matrix
    quartiles.BP <- apply(predi.mat.BP, 2, function(x) quantile(x,quant)) %>% t() 
    quartiles.BP = data.frame(quartiles.BP)
    quartiles.MFM <- apply(predi.mat.MFM, 2, function(x) quantile(x,quant)) %>% t() 
    quartiles.MFM = data.frame(quartiles.MFM)
    quartiles.DPM <- apply(predi.mat.DPM, 2, function(x) quantile(x,quant)) %>% t()
    quartiles.DPM = data.frame(quartiles.DPM)
    
    # Mudar os nomes dependendo do quantil
    colnames(quartiles.BP) = colnames(quartiles.MFM) = colnames(quartiles.DPM) = c("Q.025","Q.975")
    
    # colnames(quartiles.BP) = c("Q.GMM.025","Q.GMM.975")
    # colnames(quartiles.MFM) = c("Q.MFM.025","Q.MFM.975")
    # colnames(quartiles.DPM) = c("Q.DPM.025","Q.DPM.975")
    
    # quant.all = rbind(quartiles.BP,quartiles.MFM,quartiles.DPM)
    
    # data.BP = data.frame(y = data.distib$y.seq, true.dens = mixtures_true(data.distib$y.seq,index),quartiles.BP)
    # data.MFM = data.frame(y = data.distib$y.seq, true.dens = mixtures_true(data.distib$y.seq,index),quartiles.MFM)
    # data.DPM = data.frame(y = data.distib$y.seq, true.dens = mixtures_true(data.distib$y.seq,index),quartiles.DPM)
    
    method <- rep(c("BP", "MFM", "DPM"), each = length(y.seq))
    
    if(index != FALSE){
      data.dens <- data.frame(x = rep(y.seq,3), 
                              y = rep(true.dens,3),
                              upper = c(quartiles.BP$Q.975, quartiles.MFM$Q.975, quartiles.DPM$Q.975),
                              lower = c(quartiles.BP$Q.025, quartiles.MFM$Q.025, quartiles.DPM$Q.025),
                              method = method)
    }else{
      data.dens <- data.frame(x = rep(y.seq,3), 
                              y = quartiles.BP$Q.975,
                              upper = c(quartiles.BP$Q.975, quartiles.MFM$Q.975, quartiles.DPM$Q.975),
                              lower = c(quartiles.BP$Q.025, quartiles.MFM$Q.025, quartiles.DPM$Q.025),
                              pred  = c(mean.pred.BP,mean.pred.MFM,mean.pred.DPM),
                              method = method)
    }
    
    p1 = ggplot(data.dens, aes(x = x, y = y)) +
      geom_ribbon(aes(ymin = lower, ymax = upper, fill = method), alpha = 0.3) +
      # geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey", alpha = 0.4) +
      # geom_line() +
      # facet_wrap(~method) +
      # scale_fill_manual(values = c("Method A" = "blue", "Method B" = "green", "Method C" = "red")) +
      theme_minimal() +
      labs(title = "Three Different Ribbons for the Same Function",
           x = "X-axis",
           y = "Y-axis",
           fill = "Method") +
      theme(legend.position="bottom")
    
    if(index != FALSE){
      p1 = p1 + geom_line() +
        facet_wrap(~method) 
    }else{
      p1 = p1 + 
        geom_line(aes(x = x, y = pred, col = method)) + 
        facet_wrap(~method) 
      # geom_line(aes(x = x, y = pred), col = "red") 
    }
    
    # p2 = ggplot(data.dens, aes(x = x, y = y)) +
    #   geom_ribbon(aes(ymin = lower, ymax = upper, fill = method), alpha = 0.3) +
    #   geom_line() +
    #   facet_wrap(~method) +
    #   # scale_fill_manual(values = c("Method A" = "blue", "Method B" = "green", "Method C" = "red")) +
    #   theme_minimal() +
    #   labs(title = "Three Different Ribbons for the Same Function",
    #        x = "X-axis",
    #        y = "Y-axis",
    #        fill = "Method")
    
    p = p1
    
  }
  
  return(p)
  
}


# possum.BP = possum.BP.1
# possum.MFM = possum.MFM.1
# possum.DPM = possum.DPM.1
# 
# y.data = y.data.1
# y.data = y.data.lap

plot.possum.unc.quant.dc = function(possum.BP, possum.MFM, possum.DPM, K.sel, y.data, scale.plot = FALSE,
                                    index.possum = FALSE,index.pred = FALSE, quant = c(0.025,0.975)){
  
  param.all = rbind(possum.BP[[1]],possum.MFM[[1]],possum.DPM[[1]])
  
  #------------------------------------#
  
  if(index.possum != FALSE){
    
    mu.seq = paste0("mu", 1:K.sel)
    sigma.seq = paste0("sigma", 1:K.sel)
    prob.seq = paste0("prob", 1:K.sel)
    
    if(index.possum == "example_02"){
      data.lines.mu = data.frame(vec = c(19,19,23,29,33), variable = mu.seq, param = "mu")
      data.lines.sigma = data.frame(vec = c(5,1,1,1,2), variable = sigma.seq, param = "sigma")
      data.lines.prob = data.frame(vec = c(0.2,0.2,0.25,0.2,0.15), variable = prob.seq, param = "prob")
    }else if(index.possum == "examp_laplace"){
      data.lines.mu = data.frame(vec = c(-5,5), variable = mu.seq, param = "mu")
      data.lines.sigma = data.frame(vec = c(1.5,1), variable = sigma.seq, param = "sigma")
      data.lines.prob = data.frame(vec = c(0.4,0.6), variable = prob.seq, param = "prob")
    }
    data.lines = rbind(data.lines.mu,data.lines.sigma,data.lines.prob)
  }
  
  data_summary <- function(x){
    m <- mean(x)
    # ymin <- m-sd(x)
    # ymax <- m+sd(x)
    ymin <- as.numeric(quantile(x,probs = 0.975))
    ymax <- as.numeric(quantile(x,probs = 0.025))
    return(c(y=m,ymin=ymin,ymax=ymax))
  }
  
  p1 = ggplot(param.all) +
    # geom_violin(aes(x = method, y = value),alpha = 0.5, draw_quantiles = c(0.025,0.975)) +
    geom_violin(aes(x = method, y = value), alpha = 0.5) +
    facet_wrap(~ variable, scales = "free", ncol = K.sel) +
    labs(title = "Density Plots", x = "Value", y = "Density") +
    theme_minimal() +
    # geom_boxplot(aes(x = method, y = value),width=0.1, outlier.shape = NA) +
    # stat_summary(aes(x = method, y = value), fun.y=mean, geom="point", size=2, color="red", shape = 4) +
    # stat_summary(aes(x = method, y = value),fun.data=mean_sdl, mult=1,  geom="pointrange", color="red")
    stat_summary(aes(x = method, y = value),fun.data=data_summary, color="black", size=0.2, shape = 20) +
    # geom_hline(data = data.lines, aes(yintercept = vec), col = "red") +
    # coord_flip() +
    theme(legend.position="bottom")
  
  # p1 = p1 + geom_vline(data = data.points, aes(xintercept = vec, col = method))
  
  if(index.possum != FALSE){
    # p1 = p1 + geom_vline(data = data.lines, aes(xintercept = vec), col = "red")
    p1 = p1 + geom_hline(data = data.lines, aes(yintercept = vec), col = "red")
  }
  
  param.sum = p1
  
  #------------------------------------#
  #------------------------------------#
  
  if(sum(scale.plot) == FALSE){
    min.y<-min(y.data)-0.5*sqrt(var(y.data))
    max.y<-max(y.data)+0.5*sqrt(var(y.data))
  }else{
    min.y = scale.plot[1]
    max.y = scale.plot[2]
  }
  
  y.seq = seq(min.y, max.y, length.out = 500)
  
  # if(index.possum != FALSE){
  #   true.dens = mixtures_true(y.seq,index.possum)
  # }  
  
  if(index.possum != FALSE){
    if(index.possum == "example_01"){
      true.dens =  0.5*dnorm(y.seq,mean = 2,sd = 1) + 0.5*dnorm(y.seq,mean = 6,sd = 1)
    }else if(index.possum == "example_02"){
      true.dens =  0.2*dnorm(y.seq,mean = 19,sd = sqrt(5)) + 0.2*dnorm(y.seq,mean = 19,sd = 1) +
        0.25*dnorm(y.seq,mean = 23, sd = 1) + 0.2*dnorm(y.seq,mean = 29, sd = 1) + 0.15*dnorm(y.seq,mean = 33, sd = sqrt(2))
    }else if(index.possum == "examp_laplace"){
      true.dens = 0.4*LaplacesDemon::dlaplace(y.seq,location = -5,scale = 1.5) + 0.6*LaplacesDemon::dlaplace(y.seq,location = 5, scale  = 1)
    }
  }
  
  list.param.BP = list(mean = possum.BP[[3]]$mean_post, 
                       sigmasq = possum.BP[[3]]$sigmasq_post, 
                       pro = possum.BP[[3]]$prob_post)
  
  list.param.MFM = list(mean = possum.MFM[[3]]$mean_post, 
                        sigmasq = possum.MFM[[3]]$sigmasq_post, 
                        pro = possum.MFM[[3]]$prob_post)
  
  list.param.DPM = list(mean = possum.DPM[[3]]$mean_post, 
                        sigmasq = possum.DPM[[3]]$sigmasq_post, 
                        pro = possum.DPM[[3]]$prob_post)
  
  N = nrow(list.param.BP$mean)
  
  predi.mat.BP = matrix(0,N,length(y.seq))
  for(i in 1:N){
    for(j in 1:length(y.seq)){
      predi.mat.BP[i,j] = sum(list.param.BP$pro[i,]*dnorm(y.seq[j],mean = list.param.BP$mean[i,],sd = sqrt(list.param.BP$sigmasq[i,])))
    }
    # cat(i)
  }
  
  predi.mat.MFM = matrix(0,N,length(y.seq))
  for(i in 1:N){
    for(j in 1:length(y.seq)){
      predi.mat.MFM[i,j] = sum(list.param.MFM$pro[i,]*dnorm(y.seq[j],mean = list.param.MFM$mean[i,],sd = sqrt(list.param.MFM$sigmasq[i,])))
    }
    # cat(i)
  }
  
  predi.mat.DPM = matrix(0,N,length(y.seq))
  for(i in 1:N){
    for(j in 1:length(y.seq)){
      predi.mat.DPM[i,j] = sum(list.param.DPM$pro[i,]*dnorm(y.seq[j],mean = list.param.DPM$mean[i,],sd = sqrt(list.param.DPM$sigmasq[i,])))
    }
    # cat(i)
  }
  
  mean.pred.BP = apply(predi.mat.BP, 2, mean)
  mean.pred.MFM = apply(predi.mat.MFM, 2, mean)
  mean.pred.DPM = apply(predi.mat.DPM, 2, mean)
  
  # Apply the function to each column of the matrix
  quartiles.BP <- apply(predi.mat.BP, 2, function(x) quantile(x,quant)) %>% t() 
  quartiles.BP = data.frame(quartiles.BP)
  
  quartiles.MFM <- apply(predi.mat.MFM, 2, function(x) quantile(x,quant)) %>% t() 
  quartiles.MFM = data.frame(quartiles.MFM)
  
  quartiles.DPM <- apply(predi.mat.DPM, 2, function(x) quantile(x,quant)) %>% t()
  quartiles.DPM = data.frame(quartiles.DPM)
  
  colnames(quartiles.BP) = colnames(quartiles.MFM) = colnames(quartiles.DPM) = c("lower","upper")
  
  method <- rep(c("BP", "MFM", "DPM"), each = length(y.seq))
  
  # optim.pred.BP = mixt.dss(y.seq, list.param = possum.BP$point.estim, sel.K = possum.BP$K.sel)
  # optim.pred.MFM = mixt.dss(y.seq, list.param = possum.MFM$point.estim, sel.K = possum.MFM$K.sel)
  # optim.pred.DPM = mixt.dss(y.seq, list.param = possum.DPM$point.estim, sel.K = possum.DPM$K.sel)
  
  setwd("/Users/hbolfarine/")
  fs.MFM = read.delim("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/Julia_source/output_MFM_DPM/result_dens_MFM.txt", header = FALSE)
  fs.DPM = read.delim("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/Julia_source/output_MFM_DPM/result_dens_DPM.txt", header = FALSE)
  
  if(index.pred == "examp_laplace"){
    fs.BP = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/BP_pred_func/func_pred_exemp_laplace.csv", header = TRUE)
  }else if(index.pred == "example_02"){
    fs.BP = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/BP_pred_func/BP_pred_func.csv", header = TRUE)
  }else if(index.pred == "galaxy"){
    fs.BP = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/BP_pred_func/func_galaxy.csv", header = TRUE)
  }
  # optim.pred.BP = mixt.dss(y.seq, list.param = possum.BP$point.estim, sel.K = possum.BP$K.sel)
  # optim.pred.MFM = mixt.dss(y.seq, list.param = fs.MFM, sel.K = possum.MFM$K.sel)
  # optim.pred.DPM = mixt.dss(y.seq, list.param = fs.DPM, sel.K = possum.DPM$K.sel)
  
  optim.pred.BP = fs.BP$x
  optim.pred.MFM = fs.MFM$V1
  optim.pred.DPM = fs.DPM$V1
  
  if(index.possum != FALSE){
    data.dens <- data.frame(x = rep(y.seq,3), 
                            y = rep(true.dens,3),
                            upper = c(quartiles.BP$upper, quartiles.MFM$upper, quartiles.DPM$upper),
                            lower = c(quartiles.BP$lower, quartiles.MFM$lower, quartiles.DPM$lower),
                            pred  = c(mean.pred.BP,mean.pred.MFM,mean.pred.DPM),
                            estim = c(optim.pred.BP,optim.pred.MFM,optim.pred.DPM),
                            method = method)
  }else{
    data.dens <- data.frame(x = rep(y.seq,3), 
                            y = c(mean.pred.BP,mean.pred.MFM,mean.pred.DPM),
                            upper = c(quartiles.BP$upper, quartiles.MFM$upper, quartiles.DPM$upper),
                            lower = c(quartiles.BP$lower, quartiles.MFM$lower, quartiles.DPM$lower),
                            pred  = c(mean.pred.BP,mean.pred.MFM,mean.pred.DPM),
                            estim = c(optim.pred.BP,optim.pred.MFM,optim.pred.DPM),
                            method = method)
  }
  
  p2 = ggplot(data.dens, aes(x = x, y = y)) +
    geom_histogram(data = as.data.frame(y.data), aes(x = y.data, y = ..density..), bins = 20 ,colour="darkgrey", fill="white") +
    # geom_ribbon(aes(ymin = lower, ymax = upper, fill = method), alpha = 0.4) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
    # geom_line(aes(x = x, y = pred, col = method),linewidth = 0.3) + 
    geom_line(aes(x = x, y = pred),linewidth = 0.4, col = "red") + 
    geom_line(aes(x = x, y = estim),  linetype = "dashed", linewidth = 0.3) +
    # geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey", alpha = 0.4) +
    facet_wrap(~method) +
    # scale_fill_manual(values = c("Method A" = "blue", "Method B" = "green", "Method C" = "red")) +
    # theme_minimal() +
    theme_bw() +
    # theme_classic() +
    # theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
    labs(title = "Posterior density approximation",
         x = "Observations",
         y = "Density",
         fill = "Method") +
    theme(legend.position="bottom")
  
  if(index.possum != FALSE){
    p2 = p2 + geom_line(linewidth = 0.3)
  }
  
  # p2 = ggplot(data.dens, aes(x = x, y = y)) +
  #   geom_ribbon(aes(ymin = lower, ymax = upper, fill = method), alpha = 0.3) +
  #   geom_line() +
  #   facet_wrap(~method) +
  #   # scale_fill_manual(values = c("Method A" = "blue", "Method B" = "green", "Method C" = "red")) +
  #   theme_minimal() +
  #   labs(title = "Three Different Ribbons for the Same Function",
  #        x = "X-axis",
  #        y = "Y-axis",
  #        fill = "Method")
  
  dens.summ = p2
  if(index.possum != FALSE){
    hellinger = list(td = true.dens, m.BP = mean.pred.BP, m.MFM = mean.pred.MFM, m.DPM = mean.pred.DPM, o.BP = optim.pred.BP, o.MFM = optim.pred.MFM, o.DPM = optim.pred.DPM)
    list.plot = list(param.sum = param.sum, dens.summ = dens.summ, helling = hellinger)
  }else{
    list.plot = list(param.sum = param.sum, dens.summ = dens.summ)
  }
  
  return(list.plot)
  
}

# plot_density = "simulation"

#---------------------------------------------------------#
#---------------------------------------------------------#

# Plot density and 95% uncertainty interval for the density for univariate models 
# Simulated models 

# possum.BP[[1]]
# possum.set = possum.MFM.galaxy
# y.data = y.data.app
# K.sel = 4
# possum.set = DPM.examp_02
# y.data = data_examp_02
plot.possum.quant = function(possum.set, K.sel, y.data, scale.plot = FALSE,
                             index.possum = FALSE, index.pred = FALSE, quant = c(0.025,0.975), model = ""){
  cat("here1")
  param.all = possum.set[[1]]
  method = unique(possum.set[[1]]$method)
  
  #------------------------------------#
  
  if(index.possum != FALSE){
    
    if(index.possum == "example_02"){
      
      mu.seq = paste0("mu", 1:5)
      sigma.seq = paste0("sigma", 1:5)
      prob.seq = paste0("prob", 1:5)
      
      data.lines.mu = data.frame(vec = c(19,19,23,29,33), variable = mu.seq, param = "mu")
      data.lines.sigma = data.frame(vec = c(5,1,1,1,2), variable = sigma.seq, param = "sigma")
      data.lines.prob = data.frame(vec = c(0.2,0.2,0.25,0.2,0.15), variable = prob.seq, param = "prob")
    }else if(index.possum == "examp_laplace"){
      
      mu.seq = paste0("mu", 1:2)
      sigma.seq = paste0("sigma", 1:2)
      prob.seq = paste0("prob", 1:2)
      
      data.lines.mu = data.frame(vec = c(-5,5), variable = mu.seq, param = "mu")
      data.lines.sigma = data.frame(vec = c(1.5,1), variable = sigma.seq, param = "sigma")
      data.lines.prob = data.frame(vec = c(0.4,0.6), variable = prob.seq, param = "prob")
    }
    data.lines = rbind(data.lines.mu,data.lines.sigma,data.lines.prob)
  }
  
  data_summary <- function(x){
    m <- mean(x)
    # ymin <- m-sd(x)
    # ymax <- m+sd(x)
    ymin <- as.numeric(quantile(x,probs = 0.975))
    ymax <- as.numeric(quantile(x,probs = 0.025))
    return(c(y=m,ymin=ymin,ymax=ymax))
  }
  
  p1 = ggplot(param.all) +
    # geom_violin(aes(x = method, y = value),alpha = 0.5, draw_quantiles = c(0.025,0.975)) +
    geom_violin(aes(x = method, y = value), alpha = 0.5) +
    facet_wrap(~ variable, scales = "free", ncol = K.sel) +
    labs(title = "Density Plots", x = "Value", y = "Density") +
    theme_minimal() +
    # geom_boxplot(aes(x = method, y = value),width=0.1, outlier.shape = NA) +
    # stat_summary(aes(x = method, y = value), fun.y=mean, geom="point", size=2, color="red", shape = 4) +
    # stat_summary(aes(x = method, y = value),fun.data=mean_sdl, mult=1,  geom="pointrange", color="red")
    stat_summary(aes(x = method, y = value),fun.data=data_summary, color="black", size=0.2, shape = 20) +
    # geom_hline(data = data.lines, aes(yintercept = vec), col = "red") +
    # coord_flip() +
    theme(legend.position="bottom")

  # p1 = p1 + geom_vline(data = data.points, aes(xintercept = vec, col = method))

  if(index.possum != FALSE){
    # p1 = p1 + geom_vline(data = data.lines, aes(xintercept = vec), col = "red")
    p1 = p1 + geom_hline(data = data.lines, aes(yintercept = vec), col = "red")
  }
  # # cat("here2")
  param.sum = p1
  
  #------------------------------------#
  #------------------------------------#
  
  if(sum(scale.plot) == FALSE){
    min.y<-min(y.data)-0.5*sqrt(var(y.data))
    max.y<-max(y.data)+0.5*sqrt(var(y.data))
  }else{
    min.y = scale.plot[1]
    max.y = scale.plot[2]
  }
  
  y.seq = seq(min.y, max.y, length.out = 500)
  
  # if(index.possum != FALSE){
  #   true.dens = mixtures_true(y.seq,index.possum)
  # }  
  
  if(index.possum != FALSE){
    if(index.possum == "example_01"){
      true.dens =  0.5*dnorm(y.seq,mean = 2,sd = 1) + 0.5*dnorm(y.seq,mean = 6,sd = 1)
    }else if(index.possum == "example_02"){
      true.dens =  0.2*dnorm(y.seq,mean = 19,sd = sqrt(5)) + 0.2*dnorm(y.seq,mean = 19,sd = 1) +
        0.25*dnorm(y.seq,mean = 23, sd = 1) + 0.2*dnorm(y.seq,mean = 29, sd = 1) + 0.15*dnorm(y.seq,mean = 33, sd = sqrt(2))
    }else if(index.possum == "examp_laplace"){
      true.dens = 0.4*LaplacesDemon::dlaplace(y.seq,location = -5,scale = 1.5) + 0.6*LaplacesDemon::dlaplace(y.seq,location = 5, scale  = 1)
    }
  }
  # cat("here3")
  list.param.method = list(mean = possum.set[[3]]$mean_post, 
                           sigmasq = possum.set[[3]]$sigmasq_post, 
                           pro = possum.set[[3]]$prob_post)
  
  # list.param.MFM = list(mean = possum.MFM[[3]]$mean_post, 
  #                       sigmasq = possum.MFM[[3]]$sigmasq_post, 
  #                       pro = possum.MFM[[3]]$prob_post)
  # 
  # list.param.DPM = list(mean = possum.DPM[[3]]$mean_post, 
  #                       sigmasq = possum.DPM[[3]]$sigmasq_post, 
  #                       pro = possum.DPM[[3]]$prob_post)
  
  N = nrow(list.param.method$mean)
  
  predi.mat = matrix(0,N,length(y.seq))
  for(i in 1:N){
    for(j in 1:length(y.seq)){
      predi.mat[i,j] = sum(list.param.method$pro[i,]*dnorm(y.seq[j],mean = list.param.method$mean[i,],sd = sqrt(list.param.method$sigmasq[i,])))
    }
    # cat(i)
  }
  
  predi.point = numeric(length(y.seq))
  prob.point = possum.set$point.estim$pro[[K.sel]]
  mean.point = possum.set$point.estim$mean[[K.sel]]
  s2.point = possum.set$point.estim$sigmasq[[K.sel]]
  for(i in 1:length(y.seq)){
    predi.point[i] = sum(prob.point*dnorm(y.seq[i],mean = mean.point,sd = sqrt(s2.point)))
  }

  
  # predi.mat.MFM = matrix(0,N,length(y.seq))
  # for(i in 1:N){
  #   for(j in 1:length(y.seq)){
  #     predi.mat.MFM[i,j] = sum(list.param.MFM$pro[i,]*dnorm(y.seq[j],mean = list.param.MFM$mean[i,],sd = sqrt(list.param.MFM$sigmasq[i,])))
  #   }
  #   # cat(i)
  # }
  # 
  # predi.mat.DPM = matrix(0,N,length(y.seq))
  # for(i in 1:N){
  #   for(j in 1:length(y.seq)){
  #     predi.mat.DPM[i,j] = sum(list.param.DPM$pro[i,]*dnorm(y.seq[j],mean = list.param.DPM$mean[i,],sd = sqrt(list.param.DPM$sigmasq[i,])))
  #   }
  #   # cat(i)
  # }
  
  mean.pred = apply(predi.mat, 2, mean)
  # mean.pred.MFM = apply(predi.mat.MFM, 2, mean)
  # mean.pred.DPM = apply(predi.mat.DPM, 2, mean)
  
  # Apply the function to each column of the matrix
  quartiles <- apply(predi.mat, 2, function(x) quantile(x,quant)) %>% t() 
  quartiles = data.frame(quartiles)
  
  colnames(quartiles) = c("lower","upper")
  
  # optim.pred.BP = mixt.dss(y.seq, list.param = possum.BP$point.estim, sel.K = possum.BP$K.sel)
  # optim.pred.MFM = mixt.dss(y.seq, list.param = possum.MFM$point.estim, sel.K = possum.MFM$K.sel)
  # optim.pred.DPM = mixt.dss(y.seq, list.param = possum.DPM$point.estim, sel.K = possum.DPM$K.sel)
  if(method != "BP"){
    optim.pred = mixt.dss(y.seq, list.param = possum.set$point.estim, sel.K = possum.set$K.sel)
  }
  # setwd("/Users/hbolfarine/")
  cat("here4")
  if(method == "BP"){
    if(index.pred == "examp_laplace"){
      fs = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/BP_pred_func/func_pred_exemp_laplace.csv", header = TRUE)
    }else if(index.pred == "example_02"){
      fs = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/BP_pred_func/fun_example_02.csv", header = TRUE)
    }else if(index.pred == "galaxy"){
      fs = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/BP_pred_func/func_galaxy.csv", header = TRUE)
    }else if(index.pred == "acidity"){
      fs = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/examples_data_BP/example_acidity/func_acidity.csv", header = TRUE)
    }else if(index.pred == "enzyme"){
      fs = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/examples_data_BP/example_enzyme/func_enzyme.csv", header = TRUE)
    }else if(index.pred == "sim_data"){
      fs = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/examples_data_BP/example_mix/func_mix.csv", header = TRUE)
    }
    optim.pred = fs$x
  }
  
  # plot(y.seq,optim.pred)
  
  if(method == "MFM"){
    setwd("/Users/hbolfarine/")
    fs.MFM = read.delim("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/Julia_source/output_MFM_DPM/result_dens_MFM.txt", header = FALSE)
    optim.pred = fs.MFM$V1
  }else if(method == "DPM"){
    # setwd("/Users/hbolfarine/")
    # fs.DPM = read.delim("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/Julia_source/output_MFM_DPM/result_dens_DPM.txt", header = FALSE)
    # optim.pred = fs.DPM$V1
    fs.DPM = possum.set$pred.func.f
    optim.pred = fs.DPM
  }
  
  # optim.pred.MFM = fs.MFM$V1
  # optim.pred.DPM = fs.DPM$V1
  
  if(index.possum != FALSE){
    data.dens <- data.frame(x = y.seq, 
                            y = true.dens,
                            upper = quartiles$upper,
                            lower = quartiles$lower,
                            pred  = mean.pred,
                            estim = optim.pred,
                            true = true.dens,
                            method = method,
                            point.estm = predi.point)
    data.dens.long <- reshape(
      data.dens,
      varying   = c("pred", "estim", "point.estm", "true"),
      v.names   = "value",             # name for the gathered values
      timevar   = "line",              # name for the new grouping column
      times     = c("pred", "estim", "point.estm","true"),
      direction = "long"
    )
    p.legend = ggplot(data.dens.long, aes(x = x, y = y)) +
      geom_histogram(data = as.data.frame(y.data), aes(x = y.data, y = ..density..), bins = 30 ,colour="darkgrey", fill="white") +
      # geom_ribbon(aes(ymin = lower, ymax = upper, fill = method), alpha = 0.4) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
      geom_line(aes(y = value, color = line, linetype = line), size = 0.5) + 
      scale_color_manual(
        name   = "",
        values = c(
          pred       = "red",
          estim      = "darkgreen",
          point.estm = "blue",
          true = "black"
        ),
        labels = c(paste("Posterior expectation",model),
                   "GMM Summary estimate",
                   "Expected posterior summary",
                   "True density")
      ) +
      scale_linetype_manual(
        name   = "",
        values = c(
          pred       = "solid",
          estim      = "dashed",
          point.estm = "twodash",
          true       = "solid"
        ),
        labels = c(paste("Posterior expectation",model),
                   "GMM Summary estimate",
                   "Expected posterior summary",
                   "True density")
      ) + 
      theme_bw() +
      labs(title = "Gaussian mixture model summary",
           x = "Observations",
           y = "Density",
           fill = "Method") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.border = element_rect(color = "black", size = 0.9),
            plot.title = element_text(hjust = 0.5, size = 11),
            legend.spacing = unit(0.005, "cm"),
            legend.key.height = unit(0.2, "cm"),
            legend.title = element_text(size = 8.5),
            legend.text = element_text(size = 7.5),
            legend.background = element_rect(fill = "transparent", colour = NA),
            legend.key = element_rect(fill = "transparent", colour = NA)
            )
    
    legend = get_plot_component(p.legend, "guide-box-right", return_all = TRUE)
  }else{
    data.dens <- data.frame(x = y.seq, 
                            y = mean.pred,
                            upper = quartiles$upper,
                            lower = quartiles$lower,
                            pred  = mean.pred,
                            estim = optim.pred,
                            method = method,
                            point.estm = predi.point)
    data.dens.long <- reshape(
      data.dens,
      varying   = c("pred", "estim", "point.estm"),
      v.names   = "value",             # name for the gathered values
      timevar   = "line",              # name for the new grouping column
      times     = c("pred", "estim", "point.estm"),
      direction = "long"
    )
    p.legend = ggplot(data.dens.long, aes(x = x, y = y)) +
      geom_histogram(data = as.data.frame(y.data), aes(x = y.data, y = ..density..), bins = 30 ,colour="darkgrey", fill="white") +
      # geom_ribbon(aes(ymin = lower, ymax = upper, fill = method), alpha = 0.4) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
      geom_line(aes(y = value, color = line, linetype = line), size = 0.5) + 
      scale_color_manual(
        name   = "",
        values = c(
          pred       = "red",
          estim      = "darkgreen",
          point.estm = "blue"
        ),
        labels = c(paste("Posterior expectation",model),
                   "GMM Summary estimate",
                   "Expected posterior summary")
      ) +
      scale_linetype_manual(
        name   = "",
        values = c(
          pred       = "solid",
          estim      = "dashed",
          point.estm = "twodash"
        ),
        labels = c(paste("Posterior expectation",model),
                   "GMM Summary estimate",
                   "Expected posterior summary")
      ) + 
      theme_bw() +
      labs(title = "Gaussian mixture model summary",
           x = "Observations",
           y = "Density",
           fill = "Method") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.border = element_rect(color = "black", size = 0.9),
            plot.title = element_text(hjust = 0.5, size = 11),
            legend.spacing = unit(0.005, "cm"),
            legend.key.height = unit(0.2, "cm"),
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 7.5),
            legend.background = element_rect(fill = "transparent", colour = NA),
            legend.key = element_rect(fill = "transparent", colour = NA))
    
    legend = get_plot_component(p.legend, "guide-box-right", return_all = TRUE)
  }
  

  # data.dens = data.dens.long

  # grid::grid.draw(legend)
  
  # p2
  
  # cat("here5")
  p2 = ggplot(data.dens, aes(x = x, y = y)) +
    geom_histogram(data = as.data.frame(y.data), aes(x = y.data, y = ..density..), bins = 30 ,colour="darkgrey", fill="white") +
    # geom_ribbon(aes(ymin = lower, ymax = upper, fill = method), alpha = 0.4) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    geom_line(aes(x = x, y = point.estm),linewidth = 0.4, col = "blue", linetype = "twodash") +
    geom_line(aes(x = x, y = pred),linewidth = 0.4, col = "red") +
    geom_line(aes(x = x, y = estim), linetype = "dashed", linewidth = 0.4) +
    # geom_point(data = as.data.frame(y.data), aes(x = y.data, y = -0.02), color = "blue", shape = 16) +
    # geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey", alpha = 0.4) +
    # facet_wrap(~method) +
    # scale_fill_manual(values = c("Method A" = "blue", "Method B" = "green", "Method C" = "red")) +
    theme_bw() +
    # axis.text = element_text(size = 14),
    # axis.title.x = element_text(size = 16),
    # axis.title.y = element_text(size = 16)) +
    # labs(title = "Posterior density approximation",
    labs(title = "Gaussian mixture model summary",
         x = "Observations",
         y = "Density",
         fill = "Method") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", size = 0.9),
          plot.title = element_text(hjust = 0.5, size = 11)) #,
  
  if(index.possum != FALSE){
    p2 = p2 + geom_line(linewidth = 0.3)
  }
  
  
  dens.summ = p2
  if(index.possum != FALSE){
    hellinger = list(td = true.dens, mean.pred = mean.pred, optim.pred = optim.pred)
    list.plot = list(param.sum = param.sum, dens.summ = dens.summ, helling = hellinger, legend = legend)
  }else{
    list.plot = list(param.sum = param.sum, dens.summ = dens.summ, legend = legend)
  }
  
  return(list.plot)
  
}

#---------------------------------------------------------#
#---------------------------------------------------------#

mixtures_true_mult = function(x,y,index){
  if(index == 1){
    vec = as.numeric(c(x,y))
    dens_original = 0.5*dmvnorm(vec, mean = c(10, 10), sigma = matrix(c(1, 0, 0, 1), ncol = 2)) + 
      0.5*dmvnorm(vec, mean = c(0, 5), sigma = matrix(c(1, 0, 0, 1), ncol = 2))
    
    return(dens_original)
    
  }else if(index == "example_02"){
    vec = as.numeric(c(x,y))
    mean_vector2 = c(7, 4) 
    covariance_matrix_01 <- matrix(c(2.5, 0, 0, 0.2), ncol = 2)
    Rotate_matrix <- matrix(c(cos(pi/4),sin(pi/4),-sin(pi/4),cos(pi/4)), ncol = 2)
    covariance_matrix_02 <- Rotate_matrix%*%covariance_matrix_01%*%t(Rotate_matrix)
    mean_vector3 = c(6, 2) 
    covariance_matrix_03 <- matrix(c(3, 0, 0, 0.1), ncol = 2)
    
    dens_original = 0.45*dmvnorm(vec, mean = c(4, 4) , sigma = matrix(c(1, 0, 0, 1), ncol = 2)) + 
      0.3*dmvnorm(vec, mean = mean_vector2, sigma = covariance_matrix_02) +
      0.25*dmvnorm(vec, mean = mean_vector3, sigma = covariance_matrix_03)
    
    return(dens_original)
  }
}

#---------------------------------------------------------#

mixt.dss.mult = function(x,y,list.param,sel.K){
  vec = as.numeric(c(x,y))
  sum.final = 0
  sum.temp  <- numeric(sel.K) 
  for(i in 1:sel.K){
    # cat(i)
    sum.temp[i] = list.param$pro[[sel.K]][i]*dmvnorm(vec,mean = list.param$mean[[sel.K]][,i],sigma = list.param$sigmasq[[sel.K]][,,i])
  }
  sum.final = sum(sum.temp)
  return(sum.final)
}

#---------------------------------------------------------#

plot_pred_dens_mult = function(SFM.dens,DPM.dens,MFM.dens, y.data, index = "example_02", sel.K = 1, plot_density = FALSE){
  
  data.bind = rbind(SFM.dens[[1]],DPM.dens[[1]],MFM.dens[[1]])
  data.bind$method = factor(data.bind$method,levels = c("DPM.mult","MFM.mult","SFM.mult"))
  
  # Plot densities 
  # Old Faithful
  if(plot_density == "simulation"){
    cat("simulation")
    
    if(index == 1){
      # Simulation 1
      x.dist <- seq(-3, 15, length.out = 200)
      y.dist <- seq(-3, 20, length.out = 200)
      dat <- expand_grid(x = x.dist, y = y.dist)
      
      z1 = numeric(dim(dat)[1])
      for(i in 1:dim(dat)[1]){
        z1[i] = mixtures_true_mult(x = dat[i,1], y = dat[i,2], index = index)
      }
      
    }else if(index == "example_02"){
      cat("simulation")
      x.dist <- seq(-2, 10, length.out = 200)
      y.dist <- seq(-2, 20, length.out = 200)
      dat <- expand_grid(x = x.dist, y = y.dist)
      
      z1 = numeric(dim(dat)[1])
      for(i in 1:dim(dat)[1]){
        z1[i] = mixtures_true_mult(x = dat[i,1], y = dat[i,2], index = index)
      }
    }
    
    z2 = numeric(dim(dat)[1])
    z3 = numeric(dim(dat)[1])
    z4 = numeric(dim(dat)[1])
    for(i in 1:dim(dat)[1]){
      z2[i] = mixt.dss.mult(x = dat[i,1], y = dat[i,2],list.param = SFM.dens[[2]],sel.K)
      z3[i] = mixt.dss.mult(x = dat[i,1], y = dat[i,2],list.param = MFM.dens[[2]],sel.K)
      z4[i] = mixt.dss.mult(x = dat[i,1], y = dat[i,2],list.param = DPM.dens[[2]],sel.K)
    }  
    
    data_distrib = data.frame(dat,True = z1, SFM = z2, MFM = z3, DPM = z4)
    
    data_distrib = data_distrib %>% 
      as_tibble() %>% 
      pivot_longer(-c(1,2))
    
    data_distrib$name = factor(data_distrib$name,levels = c("True","SFM","MFM","DPM"))
    
    y.data = data.frame(X = y.data[,1], Y = y.data[,2])
    
    p = ggplot(data_distrib) + 
      geom_point(data = y.data, aes(x = X,y = Y),col = "darkgray",shape = 3, size = 1) +
      stat_contour(aes(x = x, y = y, z = value, color = name, linetype=name)) +
      # stat_contour(aes(x = x, y = y, z = value, linetype=name), col = "black") +
      # theme_minimal() +
      # scale_color_manual(values = c("True" = "black", "SFM" = "blue", "MFM" = "green", "DPM" = "red")) +
      ggtitle("True density and summary estimates") +
      xlab(TeX("$y_1$")) +
      ylab(TeX("$y_2$")) +
      theme_light() +
      theme(legend.position="bottom") +
      theme(legend.position = "bottom", legend.title = element_text(size = 10), plot.title = element_text(size = 12, hjust = 0.5))
    # theme_classic()
    
    cat("here")
    return(p)
  }else if(plot_density == "applic"){
    x.dist <- seq(min(y.data[,1])-2, max(y.data[,1]), length.out = 200)
    y.dist <- seq(min(y.data[,2])-10, max(y.data[,2]), length.out = 200)
    dat <- expand_grid(x = x.dist, y = y.dist)
    
    z2 = numeric(dim(dat)[1])
    z3 = numeric(dim(dat)[1])
    z4 = numeric(dim(dat)[1])
    for(i in 1:dim(dat)[1]){
      z2[i] = mixt.dss.mult(x = dat[i,1], y = dat[i,2],list.param = SFM.dens[[2]],sel.K)
      z3[i] = mixt.dss.mult(x = dat[i,1], y = dat[i,2],list.param = MFM.dens[[2]],sel.K)
      z4[i] = mixt.dss.mult(x = dat[i,1], y = dat[i,2],list.param = DPM.dens[[2]],sel.K)
    }  
    
    data_distrib = data.frame(dat, SFM = z2, MFM = z3, DPM = z4)
    
    data_distrib = data_distrib %>% 
      as_tibble() %>% 
      pivot_longer(-c(1,2))
    
    data_distrib$name = factor(data_distrib$name,levels = c("SFM","MFM","DPM"))
    
    y.data = data.frame(X = y.data[,1], Y = y.data[,2])
    
    p = ggplot(data_distrib) + 
      geom_point(data = y.data, aes(x = X,y = Y),col = "darkgray",shape = 3, size = 1) +
      stat_contour(aes(x = x, y = y, z = value, color = name, linetype=name)) +
      # stat_contour(aes(x = x, y = y, z = value, linetype=name), col = "black") +
      # theme_minimal() +
      # scale_color_manual(values = c("True" = "black", "SFM" = "blue", "MFM" = "green", "DPM" = "red")) +
      # ggtitle("True and approximate densities") +
      ggtitle("True density and summary estimates") +
      xlab(TeX("Eruptions$")) +
      ylab(TeX("$Waiting$")) +
      # xlab(TeX("$y_1$")) +
      # ylab(TeX("$y_2$")) +
      theme_light() +
      theme(legend.position = c(0.8, 0.8), plot.title = element_text(size = 12, hjust = 0.5), 
            axis.text.x = element_text(size = 8)) +
      theme(legend.position="bottom")
    return(p)
  }
  
}

plot.possum.uncertain.quant.dc.mult = function(possum.SFM, possum.MFM, possum.DPM, K.sel, y.data, scale.plot = FALSE,
                                               index.possum = FALSE, quant = c(0.025,0.975), examp = "simulation"){
  
  # param.all = rbind(temp.SFM.mult[[1]],temp.MFM.mult[[1]],temp.DPM.mult[[1]])
  
  #------------------------------------#
  #------------------------------------#
  #------------------------------------#
  if(examp == "simulation"){
    y1.dist <- seq(-1, 12.5, length.out = 500)
    y2.dist <- seq(-1, 12.5, length.out = 500)
    
    # fy1 - 0.45*N(4,1) + 0.3*N(7,1.35) + 0.25*N(6,3)
    # fy1 - 0.45*N(4,1) + 0.3*N(4,1.35) + 0.25*N(2,0.1)
    
    true.fy1 = 0.45*dnorm(y1.dist,4,1) + 0.3*dnorm(y1.dist,7,sqrt(1.35)) + 0.25*dnorm(y1.dist,6,sqrt(3))
    true.fy2 = 0.45*dnorm(y2.dist,4,1) + 0.3*dnorm(y2.dist,4,sqrt(1.35)) + 0.25*dnorm(y2.dist,2,sqrt(0.1))
  }else if(examp == "old_faith"){
    y1.dist <- seq(-0.5, 6, length.out = 500)
    y2.dist <- seq(35, 100, length.out = 500)
    
    true.fy1 = 0.45*dnorm(y1.dist,4,1) + 0.3*dnorm(y1.dist,7,sqrt(1.35)) + 0.25*dnorm(y1.dist,6,sqrt(3))
    true.fy2 = 0.45*dnorm(y2.dist,4,1) + 0.3*dnorm(y2.dist,4,sqrt(1.35)) + 0.25*dnorm(y2.dist,2,sqrt(0.1))
  }
  # plot(y1.dist,fy1)
  # plot(y2.dist,fy2)
  
  list.param.SFM = list(mean = possum.SFM[[2]],
                        pro = possum.SFM[[1]],
                        sigmasq = possum.SFM[[3]])
  
  list.param.MFM = list(mean = possum.MFM[[2]], 
                        pro = possum.MFM[[1]],
                        sigmasq = possum.MFM[[3]])
  
  list.param.DPM = list(mean = possum.DPM[[2]], 
                        pro = possum.DPM[[1]],
                        sigmasq = possum.DPM[[3]])
  
  N = nrow(list.param.SFM[[2]])
  
  # Select only the marginal coordinates 
  sigma.vec1 = sigma.vec2 = numeric(K.sel)
  predi.mat.SFM.y1 = predi.mat.SFM.y2 = matrix(0,N,500)
  for(i in 1:N){
    for(j in 1:500){
      for(k in 1:K.sel){
        sigma.vec1[k] = diag(list.param.SFM$sigmasq[[i]][,,k])[1]
        sigma.vec2[k] = diag(list.param.SFM$sigmasq[[i]][,,k])[2]
      }
      predi.mat.SFM.y1[i,j] = sum(list.param.SFM$pro[i,]*dnorm(y1.dist[j],mean = list.param.SFM$mean[1,,i],sd = sqrt(sigma.vec1)))
      predi.mat.SFM.y2[i,j] = sum(list.param.SFM$pro[i,]*dnorm(y2.dist[j],mean = list.param.SFM$mean[2,,i],sd = sqrt(sigma.vec2)))
    }
    # cat(i)
  }
  
  # plot(y1.dist,mean.pred.SFM.y1)
  
  sigma.vec1 = sigma.vec2 = numeric(K.sel)
  predi.mat.MFM.y1 = predi.mat.MFM.y2 = matrix(0,N,500)
  for(i in 1:N){
    for(j in 1:500){
      for(k in 1:K.sel){
        sigma.vec1[k] = diag(list.param.MFM$sigmasq[[i]][,,k])[1]
        sigma.vec2[k] = diag(list.param.MFM$sigmasq[[i]][,,k])[2]
      }
      predi.mat.MFM.y1[i,j] = sum(list.param.MFM$pro[i,]*dnorm(y1.dist[j],mean = list.param.MFM$mean[1,,i],sd = sqrt(sigma.vec1)))
      predi.mat.MFM.y2[i,j] = sum(list.param.MFM$pro[i,]*dnorm(y2.dist[j],mean = list.param.MFM$mean[2,,i],sd = sqrt(sigma.vec2)))
    }
    # cat(i)
  }
  
  sigma.vec1 = sigma.vec2 = numeric(K.sel)
  predi.mat.DPM.y1 = predi.mat.DPM.y2 = matrix(0,N,500)
  for(i in 1:N){
    for(j in 1:500){
      for(k in 1:K.sel){
        sigma.vec1[k] = diag(list.param.DPM$sigmasq[[i]][,,k])[1]
        sigma.vec2[k] = diag(list.param.DPM$sigmasq[[i]][,,k])[2]
      }
      predi.mat.DPM.y1[i,j] = sum(list.param.DPM$pro[i,]*dnorm(y1.dist[j],mean = list.param.DPM$mean[1,,i],sd = sqrt(sigma.vec1)))
      predi.mat.DPM.y2[i,j] = sum(list.param.DPM$pro[i,]*dnorm(y2.dist[j],mean = list.param.DPM$mean[2,,i],sd = sqrt(sigma.vec2)))
    }
    # cat(i)
  }
  
  mean.pred.SFM.y1 = apply(predi.mat.SFM.y1, 2, mean)
  mean.pred.SFM.y2 = apply(predi.mat.SFM.y2, 2, mean)
  
  quartiles.SFM.y1 <- apply(predi.mat.SFM.y1, 2, function(x) quantile(x,quant)) %>% t() 
  colnames(quartiles.SFM.y1) = c("lower","upper")
  quartiles.SFM.y1 = data.frame(obs = y1.dist, func = true.fy1, pred = mean.pred.SFM.y1, 
                                quartiles.SFM.y1, method = "SFM", coord = "y1")
  
  quartiles.SFM.y2 <- apply(predi.mat.SFM.y2, 2, function(x) quantile(x,quant)) %>% t() 
  colnames(quartiles.SFM.y2) = c("lower","upper")
  quartiles.SFM.y2 = data.frame(obs = y2.dist, func = true.fy2, pred = mean.pred.SFM.y2, 
                                quartiles.SFM.y2, method = "SFM", coord = "y2")
  
  #------------------------------#
  
  mean.pred.MFM.y1 = apply(predi.mat.MFM.y1, 2, mean)
  mean.pred.MFM.y2 = apply(predi.mat.MFM.y2, 2, mean)
  
  quartiles.MFM.y1 <- apply(predi.mat.MFM.y1, 2, function(x) quantile(x,quant)) %>% t() 
  colnames(quartiles.MFM.y1) = c("lower","upper")
  quartiles.MFM.y1 = data.frame(obs = y1.dist, func = true.fy1, pred = mean.pred.MFM.y1, 
                                quartiles.MFM.y1, method = "MFM", coord = "y1")
  
  quartiles.MFM.y2 <- apply(predi.mat.MFM.y2, 2, function(x) quantile(x,quant)) %>% t() 
  colnames(quartiles.MFM.y2) = c("lower","upper")
  quartiles.MFM.y2 = data.frame(obs = y2.dist, func = true.fy2, pred = mean.pred.MFM.y2, 
                                quartiles.MFM.y2, method = "MFM", coord = "y2")
  
  #------------------------------#
  
  mean.pred.DPM.y1 = apply(predi.mat.DPM.y1, 2, mean)
  mean.pred.DPM.y2 = apply(predi.mat.DPM.y2, 2, mean)
  
  quartiles.DPM.y1 <- apply(predi.mat.DPM.y1, 2, function(x) quantile(x,quant)) %>% t() 
  colnames(quartiles.DPM.y1) = c("lower","upper")
  quartiles.DPM.y1 = data.frame(obs = y1.dist, func = true.fy1, pred = mean.pred.DPM.y1, 
                                quartiles.DPM.y1, method = "DPM", coord = "y1")
  
  quartiles.DPM.y2 <- apply(predi.mat.DPM.y2, 2, function(x) quantile(x,quant)) %>% t()
  colnames(quartiles.DPM.y2) = c("lower","upper")
  quartiles.DPM.y2 = data.frame(obs = y2.dist, func = true.fy2, pred = mean.pred.DPM.y2, 
                                quartiles.DPM.y2, method = "DPM", coord = "y2")
  
  #------------------------------#
  
  
  sigma.vec1 = sigma.vec2 = numeric(K.sel)
  optim.mat.SFM.y1 = optim.mat.SFM.y2 = numeric(500)
  optim.mat.MFM.y1 = optim.mat.MFM.y2 = numeric(500)
  optim.mat.DPM.y1 = optim.mat.DPM.y2 = numeric(500)
  
  for(j in 1:500){
    for(k in 1:K.sel){
      sigma.vec1[k] = diag(possum.SFM$optim.param$sigmasq[[K.sel]][,,k])[1]
      sigma.vec2[k] = diag(possum.SFM$optim.param$sigmasq[[K.sel]][,,k])[2]
    }
    optim.mat.SFM.y1[j] = sum(possum.SFM$optim.param$pro[[K.sel]]*dnorm(y1.dist[j],mean = possum.SFM$optim.param$mean[[K.sel]][1,],sd = sqrt(sigma.vec1)))
    optim.mat.SFM.y2[j] = sum(possum.SFM$optim.param$pro[[K.sel]]*dnorm(y2.dist[j],mean = possum.SFM$optim.param$mean[[K.sel]][2,],sd = sqrt(sigma.vec2)))
  }
  
  for(j in 1:500){
    for(k in 1:K.sel){
      sigma.vec1[k] = diag(possum.MFM$optim.param$sigmasq[[K.sel]][,,k])[1]
      sigma.vec2[k] = diag(possum.MFM$optim.param$sigmasq[[K.sel]][,,k])[2]
    }
    optim.mat.MFM.y1[j] = sum(possum.MFM$optim.param$pro[[K.sel]]*dnorm(y1.dist[j],mean = possum.MFM$optim.param$mean[[K.sel]][1,],sd = sqrt(sigma.vec1)))
    optim.mat.MFM.y2[j] = sum(possum.MFM$optim.param$pro[[K.sel]]*dnorm(y2.dist[j],mean = possum.MFM$optim.param$mean[[K.sel]][2,],sd = sqrt(sigma.vec2)))
  }
  
  for(j in 1:500){
    for(k in 1:K.sel){
      sigma.vec1[k] = diag(possum.DPM$optim.param$sigmasq[[K.sel]][,,k])[1]
      sigma.vec2[k] = diag(possum.DPM$optim.param$sigmasq[[K.sel]][,,k])[2]
    }
    optim.mat.DPM.y1[j] = sum(possum.DPM$optim.param$pro[[K.sel]]*dnorm(y1.dist[j],mean = possum.DPM$optim.param$mean[[K.sel]][1,],sd = sqrt(sigma.vec1)))
    optim.mat.DPM.y2[j] = sum(possum.DPM$optim.param$pro[[K.sel]]*dnorm(y2.dist[j],mean = possum.DPM$optim.param$mean[[K.sel]][2,],sd = sqrt(sigma.vec2)))
  }
  
  
  optim.mat.SFM.y1 = possum.SFM$post.dens[,1]
  optim.mat.SFM.y2 = possum.SFM$post.dens[,2]
  
  optim.SFM.y1 = data.frame(obs = y1.dist, func = true.fy1, optim = optim.mat.SFM.y1, method = "SFM", coord = "y1")
  optim.SFM.y2 = data.frame(obs = y2.dist, func = true.fy2, optim = optim.mat.SFM.y2, method = "SFM", coord = "y2")
  
  optim.mat.MFM.y1 = possum.MFM$post.dens[,1]
  optim.mat.MFM.y2 = possum.MFM$post.dens[,2]
  
  optim.MFM.y1 = data.frame(obs = y1.dist, func = true.fy1, optim = optim.mat.MFM.y1, method = "MFM", coord = "y1")
  optim.MFM.y2 = data.frame(obs = y2.dist, func = true.fy2, optim = optim.mat.MFM.y2, method = "MFM", coord = "y2")
  
  optim.mat.DPM.y1 = possum.DPM$post.dens[,1]
  optim.mat.DPM.y2 = possum.DPM$post.dens[,2]
  
  optim.DPM.y1 = data.frame(obs = y1.dist, func = true.fy1, optim = optim.mat.DPM.y1, method = "DPM", coord = "y1")
  optim.DPM.y2 = data.frame(obs = y2.dist, func = true.fy2, optim = optim.mat.DPM.y2, method = "DPM", coord = "y2")
  
  # optim.pred.SFM = mixt.dss(y.seq, list.param = possum.SFM$point.estim, sel.K = possum.BP$K.sel)
  # optim.pred.MFM = mixt.dss(y.seq, list.param = possum.MFM$point.estim, sel.K = possum.MFM$K.sel)
  # optim.pred.DPM = mixt.dss(y.seq, list.param = possum.DPM$point.estim, sel.K = possum.DPM$K.sel)
  
  
  data.bind2 <- rbind(optim.SFM.y1,optim.SFM.y2,
                      optim.MFM.y1,optim.MFM.y2,
                      optim.DPM.y1,optim.DPM.y2)
  
  data.bind = rbind(quartiles.SFM.y1,
                    quartiles.SFM.y2,
                    quartiles.MFM.y1,
                    quartiles.MFM.y2,
                    quartiles.DPM.y1,
                    quartiles.DPM.y2)
  
  data.bind = data.bind %>% 
    mutate(optim = data.bind2$optim)
  
  p2 = ggplot(data.bind, aes(x = obs, y = func)) +
    # geom_histogram(data = as.data.frame(y.data), aes(x = y.data, y = ..density..), bins = 30 ,colour="darkgrey", fill="white") +
    # geom_ribbon(aes(ymin = lower, ymax = upper, fill = method), alpha = 0.4) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.35) +
    # geom_line(aes(x = obs, y = pred, col = method),linewidth = 0.3) + 
    geom_line(aes(x = obs, y = pred),col = "red", linewidth = 0.3) + 
    geom_line(aes(x = obs, y = optim),  linetype = "dashed", linewidth = 0.3) +
    # geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey", alpha = 0.4) +
    facet_wrap(~coord + method, scales = "free", ncol = 3) +
    # scale_fill_manual(values = c("Method A" = "blue", "Method B" = "green", "Method C" = "red")) +
    theme_minimal() +
    ggtitle("Marginal density posterior summaries") +
    xlab("Observations") +
    ylab("Density") +
    theme(legend.position="bottom") +
    theme(legend.position = c(0.8, 0.8), plot.title = element_text(size = 12, hjust = 0.5), 
          axis.text.x = element_text(size = 8))
  
  if(examp == "simulation"){
    dens.summ = p2 + geom_line(linewidth = 0.3)
  }
  dens.summ = p2
  return(dens.summ)
  
}

#----------------------------------------------##----------------------------------------------#
#----------------------------------------------##----------------------------------------------#

library(ggplot2)
library(cowplot)

# Define a function to create the plot and legends
# dens_summ = plots.possum.exemp$dens.summ
# k_star = K_star
# y.c = -0.01
# y.k = -0.02
# dat = dat.clust.DPM
create_custom_plot <- function(dat, dens_summ, k_star, y.c = -0.07, y.k = -0.15,
                               y.min.plot = -0.1, y.max.plot = 1, text_plot = "Observations"){
  
  dat = dat %>% 
    mutate(optim.clust.k = optim.clust.k+k_star)
  
  # Create the main plot
  p.temp <- dens_summ +
    geom_point(data = dat, aes(x = dat, y = y.c, color = clust.prob, shape = as.factor(optim.clust)), size = 2) +
    geom_point(data = dat, aes(x = dat, y = y.k, color = clust.prob.km, shape = as.factor(optim.clust.k)), size = 2) +
    scale_shape_manual(values = c(1:(2*k_star)),
                       name = "Clusters from loss:",
                       labels = c(paste("conditional:",1:max(dat$optim.clust)),paste("k-means:",1:max(dat$optim.clust)))) +
    scale_color_distiller( name = "Probabilities", palette = "Spectral", direction = -1, limits = c(0, 0.7)) +
    guides(
      shape = guide_legend(order = 1),   # Shape legend first
      color = guide_colorbar(order = 2)  # Color legend second
    ) +
    ylim(y.min.plot,y.max.plot) +
    xlab(text_plot) + 
    theme(
      # legend.position = c(0.9, 0.6),               # Inside kstar - plot - c(0.9, 0.6)
      # legend.justification = c("right", "top"), # ,   # Anchor point
      # legend.background = element_rect(
      #   fill = alpha("white", 0.0)
      #   # color = "black",
      #   # size = 0.5))
      # ),
      legend.spacing = unit(0.005, "cm"),
      legend.key.height = unit(0.2, "cm"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7.5),
      legend.position=c(.78,.63),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.key = element_rect(fill = "transparent", colour = NA))
  # legend.position=c(.60,.55))
  
  # p.legend <- ggdraw() +
  #   draw_plot(p.temp) +  # Your main plot
  #   draw_plot(legend, x = 0.5, y = 0.7, width = 0.25, height = 0.25)
  
  #)
  # theme(legend.position = "bottom", 
  #       plot.margin = margin(5.5, 5.5, 5.5, 5.5)) +
  # labs(color = "probability") +
  # guides(shape = "none")
  
  # # Create a separate legend for Group A
  # legend_A <- ggplot(dat, aes(x = dat, y = -0.005)) +
  #   geom_point(aes(shape = as.factor(optim.clust)), size = 2) +
  #   scale_shape_manual(
  #     name = "Cond",
  #     values = c(1:k_star),
  #     labels = c(paste0(1:k_star))
  #   ) +
  #   theme_void() +
  #   theme(legend.position = "right", 
  #         legend.text = element_text(size = 7),
  #         legend.title = element_text(size = 9))
  # 
  # # Create a separate legend for Group B
  # legend_B <- ggplot(dat, aes(x = dat, y = -0.005)) +
  #   geom_point(aes(shape = as.factor(optim.clust.k)), size = 2) +
  #   scale_shape_manual(
  #     name = "Kmeans",
  #     values = c((k_star+1):(2*k_star)),
  #     labels = c(paste0(1:k_star))
  #   ) +
  #   theme_void() +
  #   theme(legend.position = "right",
  #         legend.text = element_text(size = 7),
  #         legend.title = element_text(size = 9))
  # 
  # # Extract legends
  # legend_A_grob <- get_legend(legend_A)
  # legend_B_grob <- get_legend(legend_B)
  # 
  # # Combine plot and legends
  # final_plot <- plot_grid(
  #   p.temp, 
  #   plot_grid(legend_A_grob, legend_B_grob, ncol = 2, rel_heights = c(0.5,1),
  #             hjust = 1,
  #             vjust = 1),
  #   ncol = 2,
  #   rel_widths = c(4, 1)
  # )
  
  return(p.temp)
}

# Example of how to call the function
# create_custom_plot(dat, plots.possum.exemp$dens.summ, clust.prob, clust.prob.km, optim.clust, optim.clust.k, k.star)

#-------------------#
# plot of the posterior number of components for one plot

# Define the function

# counts_table = BP.galaxy[[3]]
plot_posterior_components <- function(counts_table, K.true = FALSE, legend.title = "Posterior on the Number of Components"){
  # Check if the input is a table or a named vector
  
  groups_present <- as.numeric(names(counts_table))
  counts <- as.numeric(counts_table)
  
  
  # Create a data frame with all groups from 1 to the maximum group
  all_groups <- data.frame(Group = 1:max(groups_present))
  
  # Create a data frame for present groups and their counts
  present_df <- data.frame(Group = groups_present, Count = counts)
  
  # Merge the two data frames to include all groups, filling missing counts with 0
  df <- all_groups %>%
    left_join(present_df, by = "Group") %>%
    mutate(Count = ifelse(is.na(Count), 0, Count))
  
  # Calculate the percentage for each group
  total_count <- sum(df$Count)
  
  df <- df %>%
    mutate(Percentage = (Count / total_count) * 100)
  
  # Create the ggplot
  p <- ggplot(df, aes(x = Group, y = Percentage)) +
    # geom_bar(stat = "identity", width = 0.3, fill = "black") +  # Thin bars
    geom_point(size = 2) +                      # Points
    geom_line() +                        # Connecting line
    # geom_text(aes(label = ifelse(Percentage > 0, round(Percentage, 1), "")), 
    #           vjust = -0.5, size = 3) +                                # Labels
    theme_classic() +
    labs(
      title = legend.title,
      x = "Number of Components",
      y = "Percentage (Counts)"
    ) +
    # scale_x_continuous(breaks = 1:max(df$Group)) +                    # Ensure all groups are labeled
    theme(
      # axis.text.x = element_text(vjust = 0.5, size = 9),  # Rotate x-axis labels
      # # axis.text.y = element_text(size = 9),
      # # axis.title = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 11),
      plot.margin = unit(c(.2,.2,.2,.5), "cm") # margin(t, r, b, l )
    )
  
  if(K.true != FALSE){
    p = p + geom_vline(xintercept = K.true, colour="gray"  , linetype = "dashed")
  }
  
  # Return the plot
  return(p)
}

# Discrepancy plot - individual one plot
#############################################
# out.mult = temp.MFM.mult
plot.estim.comp.mult.indv = function(out.mult, y.data, index = 1, sel.K = 1, plot_density = FALSE, title = "", kmax = 10){
  
  data.bind = out.mult[[1]]
  
  # data.bind  = data.bind %>% 
  #   mutate(method = fct_recode(method, "DPM" = "DPM.mult", "MFM" = "MFM.mult", "SFM" = "SFM.mult"))
  
  # p = ggplot(data.bind, aes(x = num.fact, y = avg.possum, color = method)) +
  p = ggplot(data.bind, aes(x = num.fact, y = avg.possum)) +
    # p = ggplot(data.bind, aes(x = num.fact, y = avg.possum)) + 
    # geom_hline(yintercept = -0.05, color = "gray") +
    # geom_hline(yintercept = 0.05, color = "gray") +
    geom_vline(xintercept = sel.K, colour="gray", linetype = "dashed") +
    geom_hline(yintercept = 0, colour="gray") +
    scale_x_continuous(breaks = 1:kmax) +
    # geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.2, position = position_dodge(0.6)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3) +
    geom_point(aes(x = num.fact, y = avg.possum),size = 1.5, shape = 1) +
    # theme(legend.position="bottom") +) +
    theme(legend.position="bottom") +
    theme(legend.position="bottom", legend.title=element_blank()) +
    # ggtitle(title) +
    ggtitle("Posterior summary discrepancy function") +
    xlab(TeX("Surrogate dimension $k$")) +
    ylab(TeX("$d_{n}^{k}$")) +
    # geom_line(aes(x = num.fact, y = avg.possum, linetype = method)) +
    geom_line(aes(x = num.fact, y = avg.possum)) +
    # scale_colour_manual(values=c("red","green","blue")) +
    scale_x_continuous(breaks = 1:kmax) +
    scale_linetype_manual(values=c(2,3,4)) +
    scale_shape(solid = FALSE) +
    # theme_light() +
    # theme(legend.position="bottom", legend.title=element_blank()) +
    # theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.position="bottom", legend.title=element_blank(),
          axis.text.x = element_text(size = 10),  # Rotate x-axis labels
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 12, hjust = 0.5), 
          plot.margin = unit(c(.2,.2,.2,.5), "cm")) #+ # margin(t, r, b, l )
  # theme_classic() +
  # theme_light() +
  # theme_minimal() +
  # theme(legend.position="bottom") +
  # facet_wrap(~method)
  
  return(p)
}
