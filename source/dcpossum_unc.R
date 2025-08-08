library(tidyverse)
library(mclust)
library(reshape)
# 
# calculate_quartiles = function(column,quant){
#   quantiles <- quantile(column, quant)
#   return(quantiles)
# }

model_names <- c(
  "EVV",  # ellipsoidal, equal volume (*)
  "VEV",  # ellipsoidal, equal shape
  "EEV",  # ellipsoidal, equal volume and equal shape
  "VVE",  # ellipsoidal, equal orientation (*)
  "EVE",  # ellipsoidal, equal volume and orientation (*)
  "VEE",  # ellipsoidal, equal shape and orientation (*)
  "EEE",  # ellipsoidal, equal volume, shape, and orientation
  "VVI",  # diagonal, varying volume and shape
  "EVI",  # diagonal, equal volume, varying shape
  "VEI",  # diagonal, varying volume, equal shape
  "EEI",  # diagonal, equal volume and shape
  "VII",  # spherical, unequal volume
  "EII"  # spherical, equal volume
)

dlaplace <- function(x, mu, b){
  1/(2*b) * exp(-abs(x-mu)/b)
}

# Check if the mixture are the same as in the sim data file
mixtures_true = function(x, index.func = ""){
  if(index.func == "example_02"){
    dens2 =  0.2*dnorm(x,mean = 19,sd = sqrt(5)) + 0.2*dnorm(x,mean = 19,sd = 1) +
      0.25*dnorm(x,mean = 23, sd = 1) + 0.2*dnorm(x,mean = 29, sd = 1) + 0.15*dnorm(x,mean = 33, sd = sqrt(2))
    return(dens2)
  }else if(index.func == "examp_laplace"){
    # dens3 = 0.5*dlaplace(x,mu = 0,b = 1) + 0.5*dlaplace(x,mu = 8, b  = 1)
    dens3 = 0.4*dlaplace(x,mu = -5,b = 1.5) + 0.6*dlaplace(x,mu = 5, b  = 1)
    return(dens3)
  }else if(index.func == "wand"){
    0.75*dnorm(x, mean = 0, sd = 1) + 0.25*dnorm(x, mean = 3/2, sd = 1/3)
  }
}

data_summary <- function(x){
  m <- mean(x)
  # ymin <- m-sd(x)
  # ymax <- m+sd(x)
  ymin <- as.numeric(quantile(x,probs = 0.975))
  ymax <- as.numeric(quantile(x,probs = 0.025))
  return(c(y=m,ymin=ymin,ymax=ymax))
}
# 
# possum_model = BP.1

# possum_model = BP.galaxy

# possum_model = DPM.galaxy

possum.unc.quant.values = function(possum_model, K.sel = 1){
  
  # posterior predictive sample by specific sample
  y.fit = possum_model[[4]]
  
  # Number of posterior predictive samples by each sample   
  H = ncol(y.fit)
  
  # Number of posterior samples 
  N = nrow(y.fit)
  
  # Matrix where store the posterior predictive samples by each sample
  pred.fit = matrix(0,N,H)
  
  # Matrix where we will store the posterior samples 
  # for the summary given the optimal dimension
  mean_post = sigmasq_post = prob_post = matrix(0,N,K.sel)
  
  # Generate posterior around optimal dimension and parameters 
  # y.fit = y.fit[-138,]
  # hist(y.fit[i,])
  # K.sel = 3
  for(i in 1:N){
    dens <- densityMclust(y.fit[i,], G = K.sel, plot = F, verbose = FALSE, modelNames = "V")
    # dens <- densityMclust(y.fit[i,], G = K.sel, plot = F, verbose = FALSE)
    dens1 = dens$density
    
    while(is.null(dens)){
      dens <- densityMclust(y.fit[i,], G = K.sel, plot = F, verbose = FALSE, modelNames = "E")
      dens1 = dens$density
    }
    
    # Save parameters for plot.
    # 
    # order = order(dens$parameters$mean)
    # mat.order[i,] = order
    
    prob_post[i,] = dens$parameters$pro
    mean_post[i,] = dens$parameters$mean
    sigmasq_post[i,] = dens$parameters$variance$sigmasq
    pred.fit[i,] = dens1
    # cat(i)
  }
  
  # sourceCpp("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/Rcpp_DSS_Mix/possum_uni_density_12_26_2024.cpp") # Replace with your actual file path
  # 
  # result <- computeDensityMclustArma(y_fit, K_sel)
  
  list.params.possum = list(prob_post = prob_post, mean_post = mean_post, sigmasq_post = sigmasq_post )
  
  # Generate tags for the posterior around the optimal parameters 
  mu.seq = paste0("mu", 1:K.sel)
  sigma.seq = paste0("sigma", 1:K.sel)
  prob.seq = paste0("prob", 1:K.sel)
  
  #------------------------------------#
  
  # Transform matrix to data frame for the mean
  mean_post = data.frame(mean_post)
  colnames(mean_post) = mu.seq
  
  # Possum estimate for the mean 
  mu_est_post = apply(mean_post,2, mean)
  
  mean.post = mean_post
  # Add the method's name to the plots
  # mean.post = mean_post %>% 
  #   mutate(method = method.info)
  
  #------------------------------------#
  
  # Transform matrix to data frame for sigma
  sigmasq_post = data.frame(sigmasq_post)
  colnames(sigmasq_post) = sigma.seq
  
  # Possum estimate for sigma 
  sigmasq_est_post = apply(sigmasq_post, 2, mean)
  
  # # Add the method's name to the plots
  # sigmasq.post = sigmasq_post %>% 
  #   mutate(method = method.info)
  sigmasq.post = sigmasq_post
  #------------------------------------#
  
  # Transform matrix to data frame for the probabilities 
  prob_post = data.frame(prob_post)
  colnames(prob_post) = prob.seq
  prob.post = prob_post
  
  # Possum estimate for prob 
  prob_est_post = apply(prob_post, 2, mean)
  
  # # Add the method's name to the plots
  # prob.post = prob_post %>% 
  #   mutate(method = method.info)
  
  #------------------------------------#
  #------------------------------------#
  
  mean.melt = melt(mean.post)
  sigmasq.melt = melt(sigmasq.post)
  prob.melt = melt(prob.post)
  
  #------------------------------------#
  
  mean.mealt.all = mean.melt %>% 
    mutate(param = "mu")
  
  sigma.melt.all = sigmasq.melt %>% 
    mutate(param = "sigma")
  
  prob.melt.all = prob.melt %>%
    mutate(param = "prob")
  
  param.all = rbind(mean.mealt.all,sigma.melt.all,prob.melt.all)
  
  param.all = param.all %>% 
    mutate(method = unique(possum_model[[1]]$method))

  
  # Estimates of the possum
  data.estm.mu = data.frame(vec = mu_est_post, variable = mu.seq, param = "mu")
  data.estm.sigma = data.frame(vec = sigmasq_est_post, variable = sigma.seq, param = "sigma")
  data.estm.prob = data.frame(vec = prob_est_post, variable = prob.seq, param = "prob")
  
  # Estimates possum
  data.points.estm = rbind(data.estm.mu, data.estm.sigma, data.estm.prob)
  
  data.possum = list(param.possum = param.all, 
                     possum.estm = data.points.estm, 
                     list.params.possum = list.params.possum,
                     point.estim = possum_model[[2]],
                     K.sel = K.sel, pred.func.f = possum_model[[8]])
  
  return(data.possum)
  
}

# K.sel = 3
# possum_model = temp.MFM.mult
# possum_model = temp.DPM.mult
# possum_model = temp.SFM.mult
# possum_model[[2]]$pro[K.sel]
#-----------------------------------------------------------------------------------------#

possum.unc.quant.values.mult = function(possum_model, K.sel = 1){
  
  # posterior predictive sample by specific sample
  y.fit = possum_model[[4]]
  
  # Number of posterior predictive samples by each sample   
  H = dim(y.fit)[3]
  
  # Number of posterior samples 
  N = dim(y.fit)[1]
  
  # Number of dimensions 
  p.dim = dim(y.fit)[2]
  
  # # Matrix where store the posterior predictive samples by each sample
  # pred.fit = matrix(0,N,H)
  
  # Matrix where we will store the posterior samples 
  # for the summary given the optimal dimension K.sel
  # (H,p.dim,M)
  
  mean_post = array(0, dim = c(p.dim,K.sel,N))
  prob_post = matrix(0,N,K.sel)
  sigmasq_post = list()
  
  # modelName = "VVV", use = "SVD"
  # Generate posterior around optimal dimension and parameters 
  for(i in 1:N){
    dens <- Mclust(y.fit[,,i], G = K.sel,plot = F, verbose = F)
    # dens1 = dens$density
    
    mod = 1
    while(is.null(dens)){
      dens <- Mclust(y.fit.pred[,,i], G = K.sel,plot = F, verbose = F, modelNames = model_names[mod])
      mod = mod + 1
    }
    
    # Save parameters for plot.
    mean_post[,,i] = dens$parameters$mean
    prob_post[i,] = dens$parameters$pro
    sigmasq_post[[i]] = dens$parameters$variance$sigma
    cat(i/N)
  }
  
  list.params.possum = list(prob_post = prob_post, 
                            mean_post = mean_post, 
                            sigmasq_post = sigmasq_post, 
                            optim.param = possum_model[[2]],
                            post.dens = possum_model[[8]])
  
  return(list.params.possum)
  
}

