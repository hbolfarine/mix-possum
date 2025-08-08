library(bayesm)
library(mclust)
library(tidyverse)
library(MASS)
library(mvtnorm)
library(patchwork)
# devtools::install_github("Non-Contradiction/JuliaCall")
library(JuliaCall)
library(Rcpp)
library(RcppArmadillo)
library(dirichletprocess)
library(invgamma)
library(stringr)
julia_setup(JULIA_HOME = "/Applications/Julia-1.11.app/Contents/Resources/julia/bin/")
# julia_setup(JULIA_HOME = "/Applications/Julia-1.9.app/Contents/Resources/julia/bin/")

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

# Extra functions
#----------------------------------------#

mat_numeric <- function(matrix_temp){
  numerical_matrix <- matrix(0,dim(matrix_temp)[1],2)
  for(i in 1:dim(matrix_temp)[1]){
    
    temp <- matrix_temp[i,]
    temp1 <- str_extract_all(temp,"-?\\d.*,")
    temp1 <- str_extract(temp1,"-?\\d.*\\d") %>% as.numeric()
    
    temp2 <- str_extract_all(temp,",.*-?\\d")
    temp2 <- str_extract(temp2,"-?\\d.*\\d") %>% as.numeric()
    
    numerical_matrix[i,] <- c(temp1,temp2)
  }
  return(numerical_matrix)
}

#----------------------------------------#

dens.bp = function(x,post.samp.Y.temp,temp.table,k.temp){
  dim.clust = dim(temp.table[2])[1]-1
  temp = numeric(dim.clust+1)
  for(k in 1:dim.clust){
    j.temp = temp.table$clust.k[k]
    p.temp = temp.table$prob[k]
    temp[k] = p.temp*dbeta(x,j.temp,k.temp-j.temp+1)
  }
  temp[dim.clust+1] = temp.table$prob[dim.clust+1]*dbeta(x,1,1)
  
  return(sum(temp))
}

#----------------------------------------##----------------------------------------#
# Density summarization - model with dimension d = 1
#----------------------------------------##----------------------------------------#
# Bayesian nonparametric modeling in R
# We are using the DPpackage from Alejandro Jara
# y.data = y.data.app
# data.bp = "enzyme"
dcpossum.BP = function(y.data, kmax, quant.sample = 0, BP.run = F, data.bp = "",  pred.f = FALSE, scale.plot = FALSE){
  
  setwd("/Users/hbolfarine/")
  
  if(BP.run == TRUE){
    cat("Run DPpackage code on posit!!")
    write.csv(y.data,"Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/data_bp.csv")
    break
  }
  
  # Send data to were it will be loaded to the script
  # post.samp.par = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/bp_par.csv")
  # post.samp.Y = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/post_samp_Y.csv")
  # post_pred = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/post_pred.csv")
  
  if(data.bp == "example_02"){
    post.samp.par = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/examples_data_BP/example_02/bp_par_example_02.csv")
    post.samp.Y = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/examples_data_BP/example_02/post_samp_Y_example_02.csv")
    post_pred = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/examples_data_BP/example_02/post_pred_example_02.csv")
    cat("example_02")
  }else if(data.bp == "galaxy"){
    post.samp.par = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/examples_data_BP/examp_galaxy/bp_par_galaxy.csv")
    post.samp.Y = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/examples_data_BP/examp_galaxy/post_samp_Y_galaxy.csv")
    post_pred = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/examples_data_BP/examp_galaxy/post_pred_galaxy.csv")
    cat("galaxy")
  }else if(data.bp == "mix"){
    post.samp.par = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/examples_data_BP/example_mix/bp_par_mix.csv")
    post.samp.Y = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/examples_data_BP/example_mix/post_samp_Y_mix.csv")
    post_pred = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/examples_data_BP/example_mix/post_pred_mix.csv")
    cat("mix")
  }else if(data.bp == "enzyme"){
    post.samp.par = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/examples_data_BP/example_enzyme/bp_par_enzyme.csv")
    post.samp.Y = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/examples_data_BP/example_enzyme/post_samp_Y_enzyme.csv")
    post_pred = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/examples_data_BP/example_enzyme/post_pred_enzyme.csv")
    cat("enzyme")
  }else if(data.bp == "acidity"){
    post.samp.par = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/examples_data_BP/example_acidity/bp_par_acidity.csv")
    post.samp.Y = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/examples_data_BP/example_acidity/post_samp_Y_acidity.csv")
    post_pred = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/examples_data_BP/example_acidity/post_acidity.csv")
    cat("acidity")
  }else if(data.bp == "examp_lapalce"){
    post.samp.par = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/examples_data_BP/examp_laplace/bp_par_lapalce.csv")
    post.samp.Y = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/examples_data_BP/examp_laplace/post_samp_Y_laplace.csv")
    post_pred = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/examples_data_BP/examp_laplace/post_pred_laplace.csv")
    cat("examp_lapalce")
  }else if(data.bp == "sim_250"){
    post.samp.par = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/sim_bp/bp_250/bp_par_example_250.csv")
    post.samp.Y = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/sim_bp/bp_250/post_samp_Y_example_250.csv")
    post_pred = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/sim_bp/bp_250/post_pred_example_250.csv")
    cat("sim_250")
  }else if(data.bp == "sim_500"){
    post.samp.par = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/sim_bp/bp_500/bp_par_example_500.csv")
    post.samp.Y = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/sim_bp/bp_500/post_samp_Y_example_500.csv")
    post_pred = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/sim_bp/bp_500/post_pred_example_500.csv")
    cat("sim_500")
  }else if(data.bp == "sim_1000"){
    post.samp.par = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/sim_bp/bp_1000/bp_par_example_1000.csv")
    post.samp.Y = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/sim_bp/bp_1000/post_samp_Y_example_1000.csv")
    post_pred = read.csv("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/sim_bp/bp_1000/post_pred_example_1000.csv")
    cat("sim_1000")
  }
  
  # hist(post_pred$x,breaks = 100)
  
  # hist(post_pred$x, breaks = 100)
  
  # length(post_pred$x)
  
  # Number of observations (posterior sample size)
  n.obs = length(y.data)
  # n.obs = dim(post.samp.Y)[2]
  
  # Save posterior samples for the polynomial degree
  post.k = post.samp.par[,1]
  
  numb.clust = post.samp.par[,2]
  
  # Save posterior samples for alpha
  post.alpha = post.samp.par[,3]
  
  # Save posterior samples for alpha
  M <- dim(post.samp.par)[1]
  sample_size <- M
  
  # create vector for posterior predictive sample
  #-------------------------------------------##-------------------------------------------#
  #-------------------------------------------##-------------------------------------------#
  y.pred = numeric(sample_size)
  
  item.Y = list()
  prob.Y = list()
  for(i in 1:sample_size){
    
    #-------------------------------------------#
    
    k.temp = post.k[i]
    alpha.temp = post.alpha[i]
    post.samp.Y.temp = post.samp.Y[i,]
    
    #-------------------------------------------#
    
    # try to this in rcpp
    # Which component the latent variable is located 
    # clust.k = numeric(n.obs)
    # for(q in 1:n.obs){
    #   for(j in 1:k.temp){
    #     if((j-1)/k.temp < post.samp.Y.temp[q] & post.samp.Y.temp[q] <= j/k.temp){
    #       clust.k[q] = j
    #     }
    #   }
    #   # cat(q)
    # }
    
    clust.k <- findInterval(post.samp.Y.temp, seq(0, 1, length.out = k.temp + 1), rightmost.closed = TRUE)
    
    # plot(clust.k,clust.k.2)
    
    temp.k.Y = data.frame(clust.k,Y.temp = as.numeric(post.samp.Y.temp))
    
    item.Y[[i]] = temp.k.Y
    
    #-------------------------------------------#
    # Create list of tables 
    temp.table.1 = table(clust.k) %>%
      data.frame()
    temp.table.1 = apply(temp.table.1,2,as.numeric) %>% 
      as.data.frame() %>% 
      mutate(prob = Freq/(n.obs+alpha.temp))
    
    # add_row(df, name = "Bob", age = 30) - option
    temp.table.2 <- data.frame(clust.k = as.numeric(k.temp+1), Freq = 0, prob = alpha.temp/(n.obs+alpha.temp))
    temp.table <- rbind(temp.table.1, temp.table.2)
    
    prob.Y[[i]] = temp.table
    
    #-------------------------------------------#
    
    k.clus = sample(temp.table$clust.k,1,prob = temp.table$prob)
    if(k.clus > k.temp){
      y.pred[i] = rbeta(1,1,1)  
    }else{
      y.pred[i] = rbeta(1,k.clus,k.temp-k.clus+1)
      # y.pred[i] = max(temp.k.Y[which(temp.k.Y[,1] == k.clus),2])
    }
    # cat(i)
  }
  
  #-------------------------------------------##-------------------------------------------#
  #-------------------------------------------##-------------------------------------------#
  
  # hist(y.pred, breaks = 100,probability = T)
  #-------------------------------------------#
  # Apply the posterior sample to the densities 
  
  post.dens.mat = matrix(0,sample_size,sample_size)
  for(j in 1:sample_size){
    
    k.temp = post.k[j]
    alpha.temp = post.alpha[j]
    post.samp.Y.temp = as.numeric(post.samp.Y[j,])
    
    # n.obs = 200
    # clust.k = numeric(n.obs)
    # for(q in 1:n.obs){
    #   for(j2 in 1:k.temp){
    #     if((j2-1)/k.temp < post.samp.Y.temp[q] & post.samp.Y.temp[q] <= j2/k.temp){
    #       clust.k[q] = j2
    #       # obs.clust[j,q] = post.samp.Y.temp[q]
    #     }
    #   }
    # }
    
    # generate function that returns a list of tables
    clust.k <- findInterval(post.samp.Y.temp, seq(0, 1, length.out = k.temp + 1), rightmost.closed = TRUE)
    # 
    # #-------------------------------------------#
    # 
    temp.table.1 = table(clust.k) %>%
      data.frame()
    temp.table.1 = apply(temp.table.1,2,as.numeric) %>%
      as.data.frame() %>%
      mutate(prob = Freq/(n.obs+alpha.temp))
    
    # add_row(df, name = "Bob", age = 30) - option
    temp.table.2 <- data.frame(clust.k = as.numeric(k.temp+1), Freq = 0, prob = alpha.temp/(n.obs+alpha.temp))
    temp.table <- rbind(temp.table.1, temp.table.2)
    
    for(i in 1:sample_size){
      post.dens.mat[j,i] = dens.bp(y.pred[i],post.samp.Y.temp,temp.table,k.temp)
    }
    # for(i in 1:sample_size){
    #   post.dens.mat[j,i] = dens.bp(y.pred[i],post.samp.Y.temp,prob.Y[[j]],k.temp)
    # } 
    # cat(j)
  }
  
  # plot(y.pred,post.dens.mat[5,])
  
  
  #-------------------------------------------##-------------------------------------------#
  #-------------------------------------------##-------------------------------------------#
  #  Generate predictive distribution
  # Needs to be converted to data between zero and one
  # min.y = scale.plot[1]
  # max.y = scale.plot[2]
  # y.seq.orig = seq(min.y, max.y, length.out = 500)
  # 
  # left.seq<-min(y.seq.orig)
  # right.seq<-max(y.seq.orig)
  # 
  # y.seq.t <- (y.seq.orig-left.seq)/(right.seq-left.seq)
  # 
  # post.seq.dens.mat = matrix(0,M,500)
  # for(j in 1:M){
  #   
  #   k.temp = post.k[j]
  #   alpha.temp = post.alpha[j]
  #   post.samp.Y.temp = as.numeric(post.samp.Y[j,])
  #   
  #   clust.k <- findInterval(post.samp.Y.temp, seq(0, 1, length.out = k.temp + 1), rightmost.closed = TRUE)
  #   
  #   temp.table.1 = table(clust.k) %>%
  #     data.frame()
  #   temp.table.1 = apply(temp.table.1,2,as.numeric) %>%
  #     as.data.frame() %>%
  #     mutate(prob = Freq/(n.obs+alpha.temp))
  #   
  #   # add_row(df, name = "Bob", age = 30) - option
  #   temp.table.2 <- data.frame(clust.k = as.numeric(k.temp+1), Freq = 0, prob = alpha.temp/(n.obs+alpha.temp))
  #   temp.table <- rbind(temp.table.1, temp.table.2)
  #   
  #   for(i in 1:500){
  #     post.seq.dens.mat[j,i] = dens.bp(y.seq.t[i],post.samp.Y.temp,temp.table,k.temp)
  #   }
  # }
  # pred.func.f.plot = apply(post.seq.dens.mat,2,mean)*(1.0/(right.seq-left.seq))
  
  # plot(y.seq.orig,pred.func.f.plot)
  # pred.func.f.plot = apply(post.seq.dens.mat,2,mean)
  # plot(y.seq.t,pred.func.f.plot)
  
  #-------------------------------------------##-------------------------------------------#
  #-------------------------------------------##-------------------------------------------#
  
  left<-min(y.data)-0.5*sqrt(var(y.data))
  right<-max(y.data)+0.5*sqrt(var(y.data))
  # x<-(resp-left)/(right-left)
  # jacob<-1.0/(right-left)
  # grids<- left+grid*(right-left)
  
  pred.post = apply(post.dens.mat,2,mean)
  # plot(y.pred,pred.post)
  
  pred.post.domain = pred.post*(1/(right-left))
  
  # Transform y.pred [0,1] -> Reals
  y.pred.domain = left + y.pred*(right-left)
  # plot(y.pred.domain,pred.post.domain)
  
  pred.fit = matrix(0,kmax,length(y.pred.domain))
  param.save = list(prob = list(), mean = list(), sigmasq = list())
  for(i in 1:kmax){
    dens <- densityMclust(y.pred.domain, G = i, plot = F, modelNames = "V", verbose = FALSE)
    dens1 = dens$density
    
    # Save parameters for plot.
    param.save$pro[[i]] = dens$parameters$pro
    param.save$mean[[i]] = dens$parameters$mean
    param.save$sigmasq[[i]] = dens$parameters$variance$sigmasq
    pred.fit[i,] = dens1
  }
  
  # par(mfrow = c(2,5))
  # for(i in 1:kmax){
  #   plot(y.pred,pred.fit[i,], ylab = paste("k = ",i, sep = ""))
  # }
  
  # par(mfrow = c(2,5))
  # for(i in 1:kmax){
  #   plot(y.pred.domain,pred.fit[i,], ylab = paste("k = ",i, sep = ""))
  # }
  
  #-------------------------------------------------#
  
  # Discrepancy function between the predictive surrogate densities
  # and the posterior predictive desinty
  discr = matrix(0,kmax,sample_size)
  for(i in 1:kmax){
    # discr[i,] = log(pred.fit[i,])-log(pred.post)
    discr[i,] = log(pred.fit[i,])-log(pred.post.domain)
  }
  
  # hist(discr[7,])
  
  # pred.post.temp = post.dens.mat*(1/(right-left))
  # # log.pred.fit
  # log.pred = log(pred.fit)
  # discr = matrix(0,kmax,sample_size)
  # log.pred.fit = log(pred.fit)
  # for(i in 1:kmax){
  #   temp.mat = t(log.pred.fit[i,]-t(log.pred))
  #   discr[i,] = apply(temp.mat,1,mean)
  # }
  
  # Collect the 95% confidence intervals
  
  # post.dens.mat.domain = post.dens.mat*(1/(right-left))
  # 
  # temp = matrix(0,M,kmax)
  # for(i in 1:sample_size){
  #   for(k in 1:kmax){
  #     temp[i,k] = mean(log(pred.fit[k,]) - log(post.dens.mat.domain[i,]))
  #   }  
  # }
  # 
  # discr = t(temp)
  # temp1 = apply(discr,1, function(x) quantile(x,0.975))
  # temp2 = apply(discr,1, function(x) quantile(x,0.025))
  
  # temp1 = apply(discr,1, function(x) max(x[abs(x - mean(x)) < sd(x)]))
  # temp2 = apply(discr,1, function(x) min(x[abs(x - mean(x)) < sd(x)]))
  
  temp1 = apply(discr,1, function(x) mean(x) + sd(x))
  temp2 = apply(discr,1, function(x) mean(x) - sd(x))
  
  # Saving outputs from the density estimation
  #-------------------------------------------------#
  list_output <- list()
  data.possum = data.frame(upper = temp1,
                           lower = temp2,
                           avg.possum = apply(discr,1,mean),
                           method = "BP",
                           num.fact = 1:kmax)
  
  list_output[[1]] <- data.possum
  list_output[[2]] <- param.save
  list_output[[3]] <- numb.clust
  # list_output[[5]] <- pred.post 
  # list_output[[6]] <- y.pred 
  list_output[[5]] <- pred.post.domain
  list_output[[6]] <- y.pred.domain
  list_output[[7]] <- discr
  #-------------------------------------------------#
  
  # If quant.sample != 0, generate a posterior predictive sample of size H given each parameter 
  # of the posterior sample 
  
  if(quant.sample != 0){
    H = quant.sample
    y.pred.quant = matrix(0,M,H)
    for(i in 1:M){
      temp.tab = prob.Y[[i]]
      k.clus = sample(temp.tab$clust.k,H,prob = temp.tab$prob,replace = TRUE)
      k.temp = (max(temp.tab$clust.k)-1)
      for(j in 1:H){
        # k.clus = sample(temp.tab$clust.k,1,prob = temp.tab$prob)
        if(k.clus[j] > k.temp){
          y.pred.quant[i,j] = rbeta(1,1,1)  
        }else{
          # y.pred.quant.temp = rbeta(1,k.clus[j],k.temp-k.clus[j]+1)
          y.pred.quant[i,j] = rbeta(1,k.clus[j],k.temp-k.clus[j]+1)
          # y.pred.quant[i,j] = max(item.Y[[i]][which(item.Y[[i]][,1] == k.clus[j]),2])
        }
      }
    }
    list_output[[4]] <- left + y.pred.quant*(right-left)
  }
  
  # hist(y.pred.quant[3,], probability = TRUE)
  
  if(pred.f == TRUE){
    # min.y = scale.plot[1]
    # max.y = scale.plot[2]
    # 
    left<-min(y.data)-0.5*sqrt(var(y.data))
    right<-max(y.data)+0.5*sqrt(var(y.data))
    # min.y = min(y.data) - 0.1*min(y.data)
    # max.y = max(y.data) + 0.1*max(y.data)
    y.seq = seq(left, right, length.out = 500)
    # left.sq<-min(y.seq)-0.1*sqrt(var(y.seq))
    left.sq<-left
    # right.sq<-max(y.seq)+0.1*sqrt(var(y.seq))
    right.sq<-right
    y.seq.t = (y.seq-left.sq)/(right.sq-left.sq)
    post.dens.mat = matrix(0,M,500)
    for(j in 1:M){
      
      k.temp = post.k[j]
      alpha.temp = post.alpha[j]
      post.samp.Y.temp = as.numeric(post.samp.Y[j,])
      
      clust.k <- findInterval(post.samp.Y.temp, seq(0, 1, length.out = k.temp + 1), rightmost.closed = TRUE)
      
      temp.table.1 = table(clust.k) %>%
        data.frame()
      temp.table.1 = apply(temp.table.1,2,as.numeric) %>%
        as.data.frame() %>%
        mutate(prob = Freq/(n.obs+alpha.temp))
      
      # add_row(df, name = "Bob", age = 30) - option
      temp.table.2 <- data.frame(clust.k = as.numeric(k.temp+1), Freq = 0, prob = alpha.temp/(n.obs+alpha.temp))
      temp.table <- rbind(temp.table.1, temp.table.2)
      
      for(i in 1:500){
        post.dens.mat[j,i] = dens.bp(y.seq.t[i],post.samp.Y.temp,temp.table,k.temp)
      }
    }
    pred.func.f = apply(post.dens.mat,2,mean)
    # plot(y.seq.t,pred.post.domain)
    
    pred.post.domain = pred.func.f*(1/(right.sq - left.sq))
    list_output[[8]] <- pred.post.domain
    write.csv(pred.post.domain,"Dropbox/DSS_MIX/DSS_MIX_06_2024/source/BP_source/sim_bp/pred_func.csv", row.names = FALSE)
  }
  
  return(list_output)
}

#----------------------------------------##----------------------------------------#
# Density summarization - model with dimension d = 1
#----------------------------------------##----------------------------------------#
# Mixture Models With a Prior on the Number of Components - Miller & Harrison 2018 
# We are using the Julia package BayesianMixtures from Miller & Harrison 2018 for inference

dcpossum.MFM = function(y.data,kmax, quant.sample = 0, pred.f = FALSE){
  
  setwd("/Users/hbolfarine/")
  
  # Send data to were it will be loaded to the script
  write.csv(y.data,"Dropbox/DSS_MIX/DSS_MIX_06_2024/source/Julia_source/data_MFM.csv", row.names = F)
  
  # Send data file to Julia script
  # Change the script to change the options in relation to the MFM model
  julia_source("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/Julia_source/run_MFM.jl")
  
  theta_matrix <- as.matrix(read.delim("output_MFM_DPM/result_theta_MFM.txt",header = F))
  matrix_weights_master <- as.matrix(read.delim("output_MFM_DPM/result_weight_MFM.txt",header = F))
  number_of_comp <- as.matrix(read.delim("output_MFM_DPM/result_t_MFM.txt",header = F))
  
  # Posterior sample size
  M <- dim(matrix_weights_master)[2]
  
  # Sample size of the posterior predictive sample
  # I this case will be the same size as the posterior sample
  sample_size <- M
  # sample_size <- 2000
  weight_list <- vector("list", sample_size)
  for(i in 1:M){
    weight_list[[i]] = which(matrix_weights_master[,1] != 0)
  }
  sample_size_M = sum(matrix_weights_master[which(matrix_weights_master[,1] != 0),1])
  
  # Re-weight matrix weights
  matrix_weights <- matrix_weights_master/sample_size_M
  
  # Parameters list
  param_list <- vector("list", M)
  
  # Collecting the posterior sample parameters 
  for(i in 1:M){
    
    # Length generated by the mixture
    length.MFM <- sum(matrix_weights[,i] != 0)
    
    # which components are non empty
    c = which(matrix_weights_master[,1] != 0)
    
    # Theta matrix
    matrix_temp <- matrix(theta_matrix[which(matrix_weights_master[,i] != 0),i],length.MFM,1)
    theta_mat <- mat_numeric(matrix_temp)
    
    # Vector of weights
    vector_weights <- matrix_weights[which(matrix_weights_master[,i] != 0),i]
    param_list[[i]] <- list(theta_mat = theta_mat, vector_weights = vector_weights)
  }
  
  # Generate posterior predictive sample
  y.pred = numeric(sample_size)
  weight_list <- list()
  for(i in 1:sample_size){
    weigths = param_list[[i]]$vector_weights
    weight_list[[i]] <- weigths
    G <- length(weigths)
    comp = sample(c(1:G),1,prob = weigths)
    
    mat.mu <- param_list[[i]]$theta_mat[comp,1]
    mat.sigma <- param_list[[i]]$theta_mat[comp,2]
    
    y.pred[i] = rnorm(1,mean = mat.mu, sd = mat.sigma)
  }
  # hist(y.pred)
  # Generate predictive density given the samples
  # and posterior parameters 
  expec.pred.post = matrix(0,M,sample_size)
  for(j in 1:M){
    for(i in 1:sample_size){
      mat_mu <- param_list[[j]]$theta_mat[,1]
      mat.sigma <- param_list[[j]]$theta_mat[,2]
      expec.pred.post[j,i]= sum(weight_list[[j]]*dnorm(y.pred[i],mean = mat_mu,sd = mat.sigma))
    }
  }
  
  # Predictive posterior estimate
  pred.post = apply(expec.pred.post,2,mean)
  
  # Generate posterior summaries from 1 to kmax components
  pred.fit = matrix(0,kmax,sample_size)
  param.save = list(prob = list(), mean = list(), sigmasq = list())
  for(i in 1:kmax){
    dens <- densityMclust(y.pred, G = i, plot = F, modelNames = "V", verbose = FALSE)
    dens1 = dens$density
    
    # Save parameters for plot.
    param.save$pro[[i]] = dens$parameters$pro
    param.save$mean[[i]] = dens$parameters$mean
    param.save$sigmasq[[i]] = dens$parameters$variance$sigmasq
    pred.fit[i,] = dens1
  }
  
  # Discrepancy function between the predictive surrogate densities 
  # and the posterior predictive desinty
  discr = matrix(0,kmax,sample_size)
  for(i in 1:kmax){
    discr[i,] = log(pred.fit[i,])-log(pred.post)
  }
  
  # mean(log.pred[1,] - log(pred.fit[1,])) 
  
  # temp = log.pred - log(pred.fit[1,])
  
  # temp = t(t(log.pred) - log(pred.fit[10,]))
  # 
  # boxplot(apply(temp,1,mean))
  # mean(temp[1,])
  
  # log.pred.mat = log(expec.pred.post)
  # discr = matrix(0,kmax,sample_size)
  # log.pred.fit = log(pred.fit)
  # for(i in 1:kmax){
  #   # temp.mat = t(log.pred.fit[i,] - t(log.pred))
  #   temp.mat = -t(sweep(t(log.pred.mat), 1, log.pred.fit[i,]))
  #   discr[i,] = apply(temp.mat,2,mean)
  # }
  
  # hellinger
  # log.pred.mat = sqrt(log(expec.pred.post))
  # discr = matrix(0,kmax,sample_size)
  # log.pred.fit = sqrt(log(pred.fit))
  # log.pred.fit = sum(log.pred.fit)
  
  # for(i in 1:kmax){
  #   # temp.mat = t(log.pred.fit[i,] - t(log.pred))
  #   temp.mat = t(sweep(t(log.pred.mat), 1, log.pred.fit[i,]))
  #   temp.mat2 = apply(temp.mat, 1, function(x) x^2)
  #   discr[i,] = apply(temp.mat2,2,mean)
  # }
  
  # boxplot(t(discr))
  
  
  # log(expec.pred.post)
  
  # Collect the 95% confidence intervals 
  
  # temp1 = apply(discr,1, function(x) quantile(x,0.975))
  # temp2 = apply(discr,1, function(x) quantile(x,0.025))
  
  # temp1 = apply(discr,1, function(x) max(x[abs(x - mean(x)) < sd(x)]))
  # temp2 = apply(discr,1, function(x) min(x[abs(x - mean(x)) < sd(x)]))
  
  temp1 = apply(discr,1, function(x) mean(x) + sd(x))
  temp2 = apply(discr,1, function(x) mean(x) - sd(x))
  
  # Saving outputs from the density estimation
  #-------------------------------------------------#
  list_output <- list()
  data.possum = data.frame(upper = temp1,
                           lower = temp2,
                           avg.possum = apply(discr,1,mean),
                           method = "MFM",
                           num.fact = 1:kmax)
  
  list_output[[1]] <- data.possum
  list_output[[2]] <- param.save
  list_output[[3]] <- number_of_comp
  list_output[[5]] <- pred.post
  list_output[[6]] <- y.pred 
  list_output[[7]] <- discr
  
  #-------------------------------------------------#
  
  
  # If quant.sample != 0, generate a posterior predictive sample of size H given each parameter 
  # of the posterior sample 
  
  if(quant.sample != 0){
    H = quant.sample
    y.pred.quant = matrix(0,M,H)
    weight_list <- list()
    for(i in 1:M){
      weigths = param_list[[i]]$vector_weights
      G <- length(weigths)
      comp = sample(c(1:G),H,prob = weigths, replace = T)
      for(j in 1:H){
        mat.mu <- param_list[[i]]$theta_mat[comp[j],1]
        mat.sigma <- param_list[[i]]$theta_mat[comp[j],2]      
        y.pred.quant[i,j] = rnorm(1,mean = mat.mu, sd = mat.sigma)
      }
    }
    list_output[[4]] <- y.pred.quant
  }
  
  if(pred.f == TRUE){
    min.y<-min(y.data)-0.5*sqrt(var(y.data))
    max.y<-max(y.data)+0.5*sqrt(var(y.data))
    y.seq = seq(min.y, max.y, length.out = 500)
    cat("Pred_F")
    expec.pred.post.f = matrix(0,M,500)
    for(j in 1:M){
      # cat(j)
      for(i in 1:500){
        mat_mu <- param_list[[j]]$theta_mat[,1]
        mat.sigma <- param_list[[j]]$theta_mat[,2]
        expec.pred.post.f[j,i]= sum(param_list[[j]]$vector_weights*dnorm(y.seq[i],mean = mat_mu,sd = mat.sigma))
      }
    }
    
    # Predictive posterior estimate
    pred.func.f = apply(expec.pred.post.f,2,mean)
    list_output[[8]] <- pred.func.f
  }
  return(list_output)
}

#----------------------------------------##----------------------------------------#
# Density summarization - model with dimension d = 1
#----------------------------------------##----------------------------------------#
# Dirichlet Mixture Models - (DPM)
# We are using the Julia package BayesianMixtures from Miller & Harrison 2018 for inference
# y.data = y.data.app
# y.data = data_examp_02
dcpossum.DPM = function(y.data, kmax, quant.sample = 0, pred.f = TRUE){
  alpha = 1
  setwd("/Users/hb23255/")
  
  # Send data to were it will be loaded to the script
  write.csv(y.data,"Dropbox/DSS_MIX/DSS_MIX_06_2024/source/Julia_source/data_DPM.csv", row.names = F)
  
  # Send data file to Julia script
  # Change the script to change the options in relation to the MFM model
  julia_source("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/Julia_source/run_DPM.jl")
  
  theta_matrix <- as.matrix(read.delim("output_MFM_DPM/result_theta_DPM.txt",header = F))
  matrix_weights_master <- as.matrix(read.delim("output_MFM_DPM/result_weight_DPM.txt",header = F))
  number_of_comp <- as.matrix(read.delim("output_MFM_DPM/result_t_DPM.txt",header = F))
  
  # Posterior sample size
  M <- dim(matrix_weights_master)[2]
  
  # Sample size of the posterior predictive sample
  # I this case will be the same size as the posterior sample
  # sample_size <- M
  sample_size <- 1000
  weight_list <- vector("list", sample_size)
  for(i in 1:M){
    weight_list[[i]] = which(matrix_weights_master[,1] != 0)
  }
  sample_size_M = sum(matrix_weights_master[which(matrix_weights_master[,1] != 0),1])
  
  # Re-weight matrix weights
  matrix_weights <- matrix_weights_master/sample_size_M
  # View(matrix_weights_master/sample_size)
  
  # Parameters list
  param_list <- vector("list", M)
  
  for(i in 1:M){
    
    # Length generated by the mixture
    length.DPM <- sum(matrix_weights[,i] != 0)
    c = which(matrix_weights_master[,1] != 0)
    
    # Theta matrix
    matrix_temp <- matrix(theta_matrix[which(matrix_weights_master[,i] != 0),i],length.DPM,1)
    theta_mat <- mat_numeric(matrix_temp)
    
    # Vector of weights
    vector_weights <- matrix_weights[which(matrix_weights_master[,i] != 0),i]
    param_list[[i]] <- list(theta_mat = theta_mat, vector_weights = vector_weights)
  }
  
  # Generate posterior predictive sample
  
  y.pred = numeric(sample_size)
  mat.mu <- mat.sigma <- numeric(sample_size)
  weight_list <- list()
  n.obs.dpm = length(y.data)
  mu.0 = (max(y.data)-min(y.data))/2
  s.0 = (max(y.data)-min(y.data))
  for(i in 1:sample_size){
    weigths = param_list[[i]]$vector_weights
    weight_list[[i]] <- weigths
    G <- length(weigths)
    
    
    weigths.temp = c(weight_list[[i]]*(n.obs.dpm/(alpha+n.obs.dpm)),alpha/(alpha+n.obs.dpm))
    
    comp = sample(c(1:(G+1)),1,prob = weigths.temp)
    
    # y.pred.temp = numeric(1000)
    # for(i in 1:1000){
    #   mu.temp = rnorm(1,mu.0,sd = s.0)
    #   b = rgamma(1,0.2,rate = 10/(s.0^2))
    #   # b = 1
    #   sigma.temp = rgamma(1, 2, rate = b)
    #   y.pred.temp[i] = rnorm(1,mu.temp, sd = sqrt(1/sigma.temp))
    # }
    
    
    if(comp == (G+1)){
      mu.temp = rnorm(1,mu.0,sd = s.0)
      b = rgamma(1,0.2,rate = 10/(s.0^2))
      # b = 1
      sigma.temp = rgamma(1, 2, rate = b)
      y.pred[i] = rnorm(1,mu.temp, sd = sqrt(1/sigma.temp))
    }else{
      mat.mu[i] <- param_list[[i]]$theta_mat[comp,1]
      mat.sigma[i] <- param_list[[i]]$theta_mat[comp,2]
      y.pred[i] = rnorm(1,mean = mat.mu[i], sd = mat.sigma[i])      
    }
    
  }
  # hist(y.pred, probability = T)
  # Generate predictive density given the samples
  # and posterior parameters 
  expec.pred.post = matrix(0,sample_size,sample_size)
  for(j in 1:sample_size){
    for(i in 1:length(y.pred)){
      mat_mu <- param_list[[j]]$theta_mat[,1]
      mat.sigma <- param_list[[j]]$theta_mat[,2]
      expec.pred.post[j,i]= sum(weight_list[[j]]*dnorm(y.pred[i],mean = mat_mu,sd = mat.sigma))
    }
  }
  
  # Predictive posterior estimate
  pred.post = apply(expec.pred.post,2,mean)
  
  # Generate posterior summaries from 1 to kmax components
  pred.fit = matrix(0,kmax,sample_size)
  param.save = list(prob = list(), mean = list(), sigmasq = list())
  for(i in 1:kmax){
    dens <- densityMclust(y.pred, G = i, plot =T, verbose = FALSE, modelNames = "V")
    # dens <- densityMclust(y.pred, G = i, plot =T, verbose = FALSE)
    dens1 = dens$density
    
    # Save parameters for plot.
    param.save$pro[[i]] = dens$parameters$pro
    param.save$mean[[i]] = dens$parameters$mean
    param.save$sigmasq[[i]] = dens$parameters$variance$sigmasq
    pred.fit[i,] = dens1
  }
  
  # Discrepancy function between the predictive surrogate densities 
  # and the posterior predictive desinty
  # discr = matrix(0,kmax,sample_size)
  # for(i in 1:kmax){
  #   discr[i,] = log(pred.fit[i,])-log(pred.post)
  # }
  
  temp = matrix(0,M,kmax)
  for(i in 1:sample_size){
    for(k in 1:kmax){
      temp[i,k] = mean(log(pred.fit[k,]) - log(expec.pred.post[i,])) #- mean(log(pred.post))
    }  
  }
  
  
  # log.pred = log(expec.pred.post)
  # discr = matrix(0,kmax,sample_size)
  # log.pred.fit = log(pred.fit)
  # for(i in 1:kmax){
  #   temp.mat = t(log.pred.fit[i,] - t(log.pred))
  #   discr[i,] = apply(temp.mat,1,mean)
  # }
  # 
  discr = t(temp)
  # boxplot(temp)
  # avg.possum = apply(discr,1,mean)
  # plot(avg.possum)
  temp1 = apply(discr,1, function(x) quantile(x,0.975))
  temp2 = apply(discr,1, function(x) quantile(x,0.025))
  
  # temp1 = apply(discr,1, function(x) max(x[abs(x - mean(x)) < sd(x)]))
  # temp2 = apply(discr,1, function(x) min(x[abs(x - mean(x)) < sd(x)]))
  
  # temp1 = apply(discr,1, function(x) mean(x) + sd(x))
  # temp2 = apply(discr,1, function(x) mean(x) - sd(x))
  
  # Saving outputs from the density estimation
  #-------------------------------------------------#
  list_output <- list()
  data.possum = data.frame(upper = temp1,
                           lower = temp2,
                           avg.possum = apply(discr,1,mean),
                           method = "DPM",
                           num.fact = 1:kmax)
  
  
  list_output[[1]] <- data.possum
  list_output[[2]] <- param.save
  list_output[[3]] <- number_of_comp
  list_output[[5]] <- pred.post
  list_output[[6]] <- y.pred
  list_output[[7]] <- discr
  #-------------------------------------------------# 
  
  # If quant.sample != 0, generate a posterior predictive sample of size H given each parameter 
  # of the posterior sample 
  
  # if(quant.sample != 0){
  #   H = quant.sample
  #   y.pred.quant = matrix(0,M,H)
  #   weight_list <- list()
  #   for(i in 1:M){
  #     weigths = param_list[[i]]$vector_weights
  #     G <- length(weigths)
  #     comp = sample(c(1:G),H,prob = weigths, replace = T)
  #     for(j in 1:H){
  #       mat.mu <- param_list[[i]]$theta_mat[comp[j],1]
  #       mat.sigma <- param_list[[i]]$theta_mat[comp[j],2]      
  #       y.pred.quant[i,j] = rnorm(1,mean = mat.mu, sd = mat.sigma)
  #     }
  #   }
  #   list_output[[4]] <- y.pred.quant
  # }
  
  if(quant.sample != 0){
    H = quant.sample
    y.pred.quant = matrix(0,M,H)
    weight_list <- list()
    for(i in 1:M){
      weigths = param_list[[i]]$vector_weights
      weight_list[[i]] <- weigths
      G <- length(weigths)
      
      weigths.temp = c(weight_list[[i]]*(n.obs.dpm/(alpha+n.obs.dpm)),alpha/(alpha+n.obs.dpm))
      
      comp = sample(c(1:(G+1)),H,prob = weigths.temp, replace = TRUE)
      
      for(j in 1:H){
        if(comp[j] == (G+1)){
          mu.temp = rnorm(1,mu.0,sd = s.0)
          b = rgamma(1,0.2,rate = 10/(s.0^2))
          sigma.temp = rgamma(1, 2, rate = b)
          y.pred.quant[i,j] = rnorm(1,mu.temp, sd = sqrt(1/sigma.temp))
        }else{
          mat.mu <- param_list[[i]]$theta_mat[comp[j],1]
          mat.sigma <- param_list[[i]]$theta_mat[comp[j],2]      
          y.pred.quant[i,j] = rnorm(1,mean = mat.mu, sd = mat.sigma)
        }
      }
      list_output[[4]] <- y.pred.quant
    }
  }
  
  if(pred.f == TRUE){
    min.y = min(y.data) - 0.1*min(y.data)
    max.y = max(y.data) + 0.1*max(y.data)
    y.seq = seq(min.y, max.y, length.out = 500)
    cat("Pred_F")
    expec.pred.post.f = matrix(0,M,500)
    for(j in 1:M){
      # cat(j)
      for(i in 1:500){
        mat_mu <- param_list[[j]]$theta_mat[,1]
        mat.sigma <- param_list[[j]]$theta_mat[,2]
        expec.pred.post.f[j,i]= sum(param_list[[j]]$vector_weights*dnorm(y.seq[i],mean = mat_mu,sd = mat.sigma))
      }
    }
    
    # Predictive posterior estimate
    pred.func.f = apply(expec.pred.post.f,2,mean)
    list_output[[8]] <- pred.func.f
  }
  
  return(list_output)
  
}

#----------------------------------------##----------------------------------------#
# Density summarization - model with dimension d >= 2
#----------------------------------------##----------------------------------------#
# Mixture Models With a Prior on the Number of Components - Miller & Harrison 2018 
# We are using the Julia package BayesianMixtures from Miller & Harrison 2018 for inference

# y.data = y.data.2
dcpossum.MFM.mult = function(y.data, kmax ,quant.sample = 1000, 
                             possum.samp = 0, f.pred = FALSE, pred.mod = "old_faith"){
  
  setwd("/Users/hbolfarine/")
  
  # Send data to were it will be loaded to the script
  write.csv(y.data,"Dropbox/DSS_MIX/DSS_MIX_06_2024/source/Julia_source/data_MFM_mult.csv", row.names = F)
  
  # Send data file to Julia script
  # Change the script to change the options in relation to the MFM model
  julia_source("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/Julia_source/run_MFM_mult.jl")
  
  # Posterior samples from MFM - mean 
  matrix_mu <- as.matrix(read.delim("output_MFM_DPM_mult/result_mult_mu.txt",header = F))
  
  # Posterior samples from MFM - sigma
  matrix_sigma <- as.matrix(read.delim("output_MFM_DPM_mult/result_mult_sigma2.txt",header = F)) 
  
  # Posterior samples from MFM - weights
  matrix_weights <- as.matrix(read.delim("output_MFM_DPM_mult/result_weight_mult.txt",header = F))
  
  # Distribution of the number of components
  number_of_clust <- read.delim("output_MFM_DPM_mult/result_t_mult.txt",header = F)
  
  # number_of_comp <- as.matrix(read.delim("output_MFM_DPM_mult/result_z_mult.txt",header = F))
  
  clust_class <- as.matrix(read.delim("output_MFM_DPM_mult/result_mult_class_MFM.txt",header = F))
  
  # Data dimension 
  p.dim = dim(y.data)[2]
  
  # Posterior sample size
  M <- dim(matrix_weights)[2]
  
  sample_size <- quant.sample
  # Posterior predictive sample
  y.pred = matrix(0,sample_size,p.dim)
  for(i in 1:sample_size){
    comp = which(matrix_weights[,i] != 0)
    weight_list = matrix_weights[comp,i]
    comp.sample = sample(comp,1,prob = weight_list)
    mu_vec = matrix_mu[comp.sample,(p.dim*i-(p.dim-1)):(p.dim*i)]
    sigma_mat = matrix(matrix_sigma[comp.sample,((p.dim^2)*i-(p.dim^2-1)):((p.dim^2)*i)],p.dim,p.dim)
    y.pred[i,] = mvrnorm(n = 1, mu = mu_vec,Sigma = sigma_mat)
  }
  
  #----------------------------------------#
  
  # Generate predictive density given the samples
  # and posterior parameters 
  # expec.pred.post = matrix(0,sample_size,sample_size)
  # for(j in 1:sample_size){
  #   comp = which(matrix_weights[,j] != 0)
  #   for(i in 1:sample_size){
  #     density.temp = numeric(length(comp))
  #     t = 1
  #     for(k in comp){
  #       mu_vec = matrix_mu[k,(p.dim*j-(p.dim-1)):(p.dim*j)]
  #       sigma_mat = matrix(matrix_sigma[k,((p.dim^2)*j-(p.dim^2-1)):((p.dim^2)*j)],p.dim,p.dim)
  #       prob = matrix_weights[k,j]
  #       density.temp[t] <- prob*mvtnorm::dmvnorm(y.pred[i,],mean = mu_vec,sigma = sigma_mat)
  #       t = t + 1
  #     }
  #     expec.pred.post[j,i] = sum(density.temp)
  #   }
  #   # cat(round(j/M,3),"\n")
  # }
  
  setwd("/Users/hbolfarine/")
  sourceCpp("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/Rcpp_DSS_Mix/pred_MFM_mult_12_23_2024.cpp")
  cat("cpp-start")
  expec.pred.post = expec_pred_post_cpp(matrix_weights,
                                        matrix_mu,
                                        matrix_sigma,
                                        y.pred,
                                        sample_size,
                                        p.dim)
  
  # Predictive posterior estimate
  pred.post = apply(expec.pred.post,2,mean)
  
  # Generate posterior summaries from 1 to kmax components
  y.pred.data <- data.frame(y.pred)
  pred.fit = matrix(0,kmax,sample_size)
  list_names <- numeric(kmax)
  param.save = list(prob = list(), mean = list(), sigmasq = list())
  for(i in 1:kmax){
    dens <- densityMclust(y.pred.data, G = i, plot = F, modelNames = "VVV")
    # dens <- densityMclust(y.pred.data, G = i, plot = F)
    
    mod = 1
    while(is.null(dens)){
      dens <- Mclust(y.pred.data, G = i,plot = F, verbose = F, modelNames = model_names[mod])
      mod = mod + 1
    }
    
    # if(length(dens$density) == 0){
    #   dens$density = density.temp
    #   dens$parameters$pro = prob.temp
    #   dens$parameters$mean = mean.temp
    #   dens$parameters$variance$sigma = sigma.temp
    #   dens$modelName = name.temp
    #   # cat("Null")
    # }
    
    # Save parameters for plot.
    param.save$pro[[i]] = dens$parameters$pro
    param.save$mean[[i]] = dens$parameters$mean
    param.save$sigmasq[[i]] = dens$parameters$variance$sigma
    pred.fit[i,] = dens$density
    list_names[i] <- dens$modelName
    
    # density.temp = dens$density
    # prob.temp = dens$parameters$pro
    # mean.temp = dens$parameters$mean
    # sigma.temp = dens$parameters$variance$sigma
    # name.temp = dens$modelName
    
  }
  
  #----------------------------------------#
  
  # Discrepancy function between the predictive surrogate densities 
  # and the posterior predictive density
  discr = matrix(0,kmax,sample_size)
  for(i in 1:kmax){
    discr[i,] = log(pred.fit[i,])-log(pred.post)
  }
  # 
  
  # temp = matrix(0,M,kmax)
  # for(i in 1:sample_size){
  #   for(k in 1:kmax){
  #     temp[i,k] = mean(log(pred.fit[k,]) - log(expec.pred.post[i,]))
  #   }  
  # }
  
  # boxplot(temp)
  # discr = t(temp)
  # temp1 = apply(discr,1, function(x) quantile(x,0.975))
  # temp2 = apply(discr,1, function(x) quantile(x,0.025))
  
  # temp1 = apply(discr,1, function(x) max(x[abs(x - mean(x)) < sd(x)]))
  # temp2 = apply(discr,1, function(x) min(x[abs(x - mean(x)) < sd(x)]))
  
  temp1 = apply(discr,1, function(x) mean(x) + sd(x))
  temp2 = apply(discr,1, function(x) mean(x) - sd(x))
  
  list_output <- list()
  data.possum = data.frame(upper = temp1,
                           lower = temp2,
                           avg.possum = apply(discr,1,mean),
                           method = "MFM.mult",
                           model_names = list_names,
                           num.fact = 1:kmax)
  
  list_output[[1]] <- data.possum
  list_output[[2]] <- param.save
  list_output[[3]] <- number_of_clust
  list_output[[5]] <- pred.post
  list_output[[6]] <- y.pred
  list_output[[7]] <- discr
  list_output[[9]] <- clust_class
  
  cat("-Post_Samp-")
  if(possum.samp != 0){
    H = possum.samp
    # y.pred.quant = matrix(0,N,H)
    y.pred.quant = array(0,dim = c(H,p.dim,M))
    for(i in 1:M){
      comp = which(matrix_weights[,i] != 0)
      weight_list = matrix_weights[comp,i]
      comp.sample = sample(comp,H,prob = weight_list,replace = T)
      
      for(j in 1:H){
        mu_vec = matrix_mu[comp.sample[j],(p.dim*i-(p.dim-1)):(p.dim*i)]
        sigma_mat = matrix(matrix_sigma[comp.sample[j],((p.dim^2)*i-(p.dim^2-1)):((p.dim^2)*i)],p.dim,p.dim)
        y.pred.quant[j,,i] = mvrnorm(1, mu = mu_vec, Sigma = sigma_mat)
      }
    }
    list_output[[4]] <- y.pred.quant
  }
  
  cat("-fpred-")
  if(f.pred == TRUE){
    
    if(pred.mod == "simulation"){
      y1.dist <- seq(-1, 12.5, length.out = 500)
      y2.dist <- seq(-1, 12.5, length.out = 500)
    }else if(pred.mod == "old_faith"){
      y1.dist <- seq(-0.5, 6, length.out = 500)
      y2.dist <- seq(35, 100, length.out = 500)
    }
    
    dens.post.1 = matrix(0,1000,500)
    dens.post.2 = matrix(0,1000,500)
    
    for(j in 1:1000){
      comp = which(matrix_weights[,j] != 0)
      density.temp.1 = numeric(length(comp))
      density.temp.2 = numeric(length(comp))
      for(i in 1:500){
        t = 1
        for(k in comp){
          mu_vec = matrix_mu[k,(p.dim*j-(p.dim-1)):(p.dim*j)]
          sigma_mat = matrix(matrix_sigma[k,((p.dim^2)*j-(p.dim^2-1)):((p.dim^2)*j)],p.dim,p.dim)
          prob = matrix_weights[k,j]
          
          density.temp.1[t] = prob*dnorm(y1.dist[i],mean = mu_vec[1],sd = sqrt(sigma_mat[1,1]))
          density.temp.2[t] = prob*dnorm(y2.dist[i],mean = mu_vec[2],sd = sqrt(sigma_mat[2,2]))
          
          t = t + 1
        }
        dens.post.1[j,i] = sum(density.temp.1)
        dens.post.2[j,i] = sum(density.temp.2)
      }
      
      expec.dens.post.1 = apply(dens.post.1,2,mean)
      expec.dens.post.2 = apply(dens.post.2,2,mean)
    }
    list_output[[8]] <- cbind(expec.dens.post.1,expec.dens.post.2)
  }
  
  return(list_output)
  
}


#----------------------------------------##----------------------------------------#
# Density summarization - model with dimension d >= 2 
#----------------------------------------##----------------------------------------#
# Dirichlet Mixture Models - (DPM)
# We are using the Julia package BayesianMixtures from Miller & Harrison 2018 for inference

# Dirichlet Process Mixture Models - (DPM)
# y.data = y.data.2
dcpossum.DPM.mult = function(y.data, kmax, quant.sample = 1000, 
                             possum.samp = 0, f.pred = FALSE, pred.mod = "simulation"){
  
  setwd("/Users/hb23255/")
  
  # Send data to were it will be loaded to the script
  write.csv(y.data,"Dropbox/DSS_MIX/DSS_MIX_06_2024/source/Julia_source/data_DPM_mult.csv", row.names = F)
  
  # Send data file to Julia script
  # Change the script to change the options in relation to the MFM model
  julia_source("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/Julia_source/run_DPM_mult.jl")
  
  # Posterior samples from DPM - mean
  matrix_mu <- read.delim("output_MFM_DPM_mult/result_mult_mu_DPM.txt",header = F) %>%
    as.matrix()
  
  # Posterior samples from DPM - sigma
  matrix_sigma <- read.delim("output_MFM_DPM_mult/result_mult_sigma2_DPM.txt",header = F) %>%
    as.matrix()
  
  # Posterior samples from DPM - weights
  matrix_weights <- read.delim("output_MFM_DPM_mult/result_weight_mult_DPM.txt",header = F) %>%
    as.matrix()
  
  # Distribution of the number of components
  number_of_comp <- read.delim("output_MFM_DPM_mult/result_t_mult_DPM.txt",header = F)
  
  clust_class <- read.delim("output_MFM_DPM_mult/result_mult_class_DPM.txt",header = F)
  
  # Data dimension 
  p.dim = dim(y.data)[2]
  
  # Posterior sample size
  M <- dim(matrix_weights)[2]
  
  sample_size <- quant.sample
  # Posterior predictive sample
  
  # Generate posterior summaries from 1 to kmax components
  y.pred = matrix(0,sample_size,p.dim)
  for(i in 1:sample_size){
    comp = which(matrix_weights[,i] != 0)
    weight_list = matrix_weights[comp,i]
    comp.sample = sample(comp,1,prob = weight_list)
    mu_vec = matrix_mu[comp.sample,(p.dim*i-(p.dim-1)):(p.dim*i)]
    sigma_mat = matrix(matrix_sigma[comp.sample,((p.dim^2)*i-(p.dim^2-1)):((p.dim^2)*i)],p.dim,p.dim)
    y.pred[i,] = mvrnorm(n = 1, mu = mu_vec,Sigma = sigma_mat)
  }
  
  # Generate posterior summaries from 1 to kmax components
  # expec.pred.post = matrix(0,M,sample_size)
  # for(j in 1:sample_size){
  #   comp = which(matrix_weights[,j] != 0)
  #   for(i in 1:sample_size){
  #     density.temp = numeric(length(comp))
  #     t = 1
  #     for(k in comp){
  #       mu_vec = matrix_mu[k,(p.dim*j-(p.dim-1)):(p.dim*j)]
  #       sigma_mat = matrix(matrix_sigma[k,((p.dim^2)*j-(p.dim^2-1)):((p.dim^2)*j)],p.dim,p.dim)
  #       prob = matrix_weights[k,j]
  #       density.temp[t] <- prob*mvtnorm::dmvnorm(y.pred[i,],mean = mu_vec,sigma = sigma_mat)
  #       t = t + 1
  #     }
  #     expec.pred.post[j,i] = sum(density.temp)
  #   }
  #   cat(round(j/M,3),"\n")
  # }
  
  setwd("/Users/hb23255/")
  sourceCpp("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/Rcpp_DSS_Mix/pred_MFM_mult_12_23_2024.cpp")
  cat("cpp-start")
  expec.pred.post = expec_pred_post_cpp(matrix_weights,
                                        matrix_mu,
                                        matrix_sigma,
                                        y.pred,
                                        sample_size,
                                        p.dim)
  
  
  # Predictive posterior estimate
  pred.post = apply(expec.pred.post,2,mean)
  
  y.pred.data <- data.frame(y.pred)
  pred.fit = matrix(0,kmax,M)
  list_names <- numeric(kmax)
  param.save = list(prob = list(), mean = list(), sigmasq = list())
  for(i in 1:kmax){
    # dens <- densityMclust(y.pred.data, G = i, plot = F)
    dens <- densityMclust(y.pred.data, G = i, plot = F, modelNames = "VVV")
    
    mod = 1
    while(is.null(dens)){
      dens <- Mclust(y.pred.data, G = K.sel,plot = F, verbose = F, modelNames = model_names[mod])
      mod = mod + 1
    }
    
    # if(length(dens$density) == 0){
    #   dens$density = density.temp
    #   dens$parameters$pro = prob.temp
    #   dens$parameters$mean = mean.temp
    #   dens$parameters$variance$sigma = sigma.temp
    #   dens$modelName = name.temp
    #   # cat("Null")
    # }
    
    # Save parameters for plot.
    param.save$pro[[i]] = dens$parameters$pro
    param.save$mean[[i]] = dens$parameters$mean
    param.save$sigmasq[[i]] = dens$parameters$variance$sigma
    pred.fit[i,] = dens$density
    list_names[i] <- dens$modelName
    
    # density.temp = dens$density
    # prob.temp = dens$parameters$pro
    # mean.temp = dens$parameters$mean
    # sigma.temp = dens$parameters$variance$sigma
    # name.temp = dens$modelName
  }
  
  # Discrepancy function between the predictive surrogate densities 
  # and the posterior predictive density
  # discr = matrix(0,kmax,M)
  # for(i in 1:kmax){
  #   discr[i,] = log(pred.fit[i,])-log(pred.post)
  # }
  
  temp = matrix(0,M,kmax)
  for(i in 1:sample_size){
    for(k in 1:kmax){
      temp[i,k] = mean(log(pred.fit[k,]) - log(expec.pred.post[i,]))
    }  
  }
  
  # boxplot(temp)
  discr = t(temp)
  temp1 = apply(discr,1, function(x) quantile(x,0.975))
  temp2 = apply(discr,1, function(x) quantile(x,0.025))
  
  # temp1 = apply(discr,1, function(x) max(x[abs(x - mean(x)) < sd(x)]))
  # temp2 = apply(discr,1, function(x) min(x[abs(x - mean(x)) < sd(x)]))
  
  # temp1 = apply(discr,1, function(x) mean(x) + sd(x))
  # temp2 = apply(discr,1, function(x) mean(x) - sd(x))
  
  list_output <- list()
  data.possum = data.frame(upper = temp1,
                           lower = temp2,
                           avg.possum = apply(discr,1,mean),
                           method = "DPM.mult",
                           model_names = list_names,
                           num.fact = 1:kmax)
  
  
  list_output[[1]] <- data.possum
  list_output[[2]] <- param.save
  list_output[[3]] <- number_of_comp
  list_output[[5]] <- pred.post
  list_output[[6]] <- y.pred
  list_output[[7]] <- discr
  list_output[[9]] <- clust_class
  
  if(possum.samp != 0){
    H = possum.samp
    # y.pred.quant = matrix(0,N,H)
    y.pred.quant = array(0,dim = c(H,p.dim,M))
    for(i in 1:M){
      comp = which(matrix_weights[,i] != 0)
      weight_list = matrix_weights[comp,i]
      comp.sample = sample(comp,H,prob = weight_list,replace = T)
      
      for(j in 1:H){
        mu_vec = matrix_mu[comp.sample[j],(p.dim*i-(p.dim-1)):(p.dim*i)]
        sigma_mat = matrix(matrix_sigma[comp.sample[j],((p.dim^2)*i-(p.dim^2-1)):((p.dim^2)*i)],p.dim,p.dim)
        y.pred.quant[j,,i] = mvrnorm(1, mu = mu_vec, Sigma = sigma_mat)
      }
    }
    list_output[[4]] <- y.pred.quant
  }
  
  if(f.pred == TRUE){
    
    if(pred.mod == "simulation"){
      y1.dist <- seq(-1, 12.5, length.out = 500)
      y2.dist <- seq(-1, 12.5, length.out = 500)
    }else if(pred.mod == "old_faith"){
      y1.dist <- seq(-0.5, 6, length.out = 500)
      y2.dist <- seq(35, 100, length.out = 500)
    }
    
    dens.post.1 = matrix(0,1000,500)
    dens.post.2 = matrix(0,1000,500)
    
    for(j in 1:1000){
      comp = which(matrix_weights[,j] != 0)
      density.temp.1 = numeric(length(comp))
      density.temp.2 = numeric(length(comp))
      for(i in 1:500){
        t = 1
        for(k in comp){
          mu_vec = matrix_mu[k,(p.dim*j-(p.dim-1)):(p.dim*j)]
          sigma_mat = matrix(matrix_sigma[k,((p.dim^2)*j-(p.dim^2-1)):((p.dim^2)*j)],p.dim,p.dim)
          prob = matrix_weights[k,j]
          
          density.temp.1[t] = prob*dnorm(y1.dist[i],mean = mu_vec[1],sd = sqrt(sigma_mat[1,1]))
          density.temp.2[t] = prob*dnorm(y2.dist[i],mean = mu_vec[2],sd = sqrt(sigma_mat[2,2]))
          
          t = t + 1
        }
        dens.post.1[j,i] = sum(density.temp.1)
        dens.post.2[j,i] = sum(density.temp.2)
      }
      
      expec.dens.post.1 = apply(dens.post.1,2,mean)
      expec.dens.post.2 = apply(dens.post.2,2,mean)
    }
    list_output[[8]] <- cbind(expec.dens.post.1,expec.dens.post.2)
  }
  
  return(list_output)
  
}

dcpossum.DPM.dir.mult = function(y.data, kmax, possum.samp = 0, f.pred = FALSE, pred.mod = "simulation"){
  
  # old, kappa = 10, nu = 5
  # old, kappa = 1, nu = 2
  # sim, kappa = 1, nu = 2
  # thyroid, kappa = 1, nu = 5
  
  kappa0 = 1
  
  # Degrees of freedom - Wishart
  nu.prior = 5
  
  # Number of observations in the data
  n.obs.dpm.mult = dim(y.data)[1]
  
  # Number of dimensions
  p.dim = dim(y.data)[2]
  
  # Prior mean 
  mu.prior = apply(y.data, 2, mean)
  
  # Covariance prior
  # T0 = diag(ncol(y.data.2))
  T0 = cov(y.data)
  
  # List of priors
  g0Priors = list(mu0 = mu.prior,
                  Lambda = T0,
                  kappa0 = kappa0,
                  nu = nu.prior) 
  
  dp.mult = DirichletProcessMvnormal(as.matrix(y.data), g0Priors = g0Priors,
                                     alphaPriors = c(2,4))
  
  cat("DPM estimate\n")
  its = 2000
  dp.mult <- Fit(dp.mult, its, progressBar = TRUE)
  dp.mult = Burn(dp.mult, 1000)
  
  sample_size = 1000
  
  cat("Generate posterior predictive sample\n")
  y.pred = matrix(0,sample_size,p.dim)
  for(i in 1:sample_size){
    
    G <- length(dp.mult$weightsChain[[i]])
    alpha.temp = dp.mult$alphaChain[i]
    weigths.temp = c(dp.mult$weightsChain[[i]]*(n.obs.dpm.mult/(alpha.temp + n.obs.dpm.mult)),alpha.temp/(alpha.temp + n.obs.dpm.mult))
    # weigths.temp = dp.mult$weightsChain[[i]]
    
    comp = sample(c(1:(G+1)) ,1 , prob = weigths.temp, replace = TRUE)
    # comp = sample(c(1:(G)) ,1 , prob = weigths.temp, replace = TRUE)
    
    if(comp == (G+1)){
      # Remember to change the base measure
      Lamb.temp = LaplacesDemon::rinvwishart(nu.prior, solve(T0))
      y.pred[i,] = mvrnorm(n = 1, mu = mu.prior,Sigma = (1/kappa0)*Lamb.temp)
    }else{
      mu_vec = dp.mult$clusterParametersChain[[i]][1]$mu[,,comp]
      sigma_mat = dp.mult$clusterParametersChain[[i]][2]$sig[,,comp]
      y.pred[i,] = mvrnorm(n = 1, mu = mu_vec, Sigma = sigma_mat)
    }
  }
  
  cat("Generate posterior predictive density\n")
  expec.pred.post = matrix(0,sample_size,sample_size)
  for(j in 1:sample_size){
    G = length(dp.mult$weightsChain[[j]])
    dens = numeric(G)
    prob = dp.mult$weightsChain[[j]]
    for(i in 1:sample_size){
      for(k in 1:G){
        mu_vec = dp.mult$clusterParametersChain[[j]][1]$mu[,,k]
        sigma_mat = dp.mult$clusterParametersChain[[j]][2]$sig[,,k]
        dens[k] = prob[k]*mvtnorm::dmvnorm(y.pred[i,], mean = mu_vec, sigma = sigma_mat)
      } 
      expec.pred.post[j,i] = sum(dens)
    }
    cat(round(j/sample_size,2),"\n")
  }
  
  # Predictive posterior estimate
  pred.post = apply(expec.pred.post,2,mean)
  
  y.pred.data <- data.frame(y.pred)
  pred.fit = matrix(0,kmax,sample_size)
  list_names <- numeric(kmax)
  param.save = list(prob = list(), mean = list(), sigmasq = list())
  for(i in 1:kmax){
    # dens <- densityMclust(y.pred.data, G = i, plot = F)
    dens <- densityMclust(y.pred.data, G = i, plot = F, modelNames = "VVV")
    
    # mod = 1
    while(is.null(dens)){
      dens <- Mclust(y.pred.data, G = i,plot = F, verbose = F, modelNames = model_names[mod])
      # mod = mod + 1
    }
    
    # if(length(dens$density) == 0){
    #   dens$density = density.temp
    #   dens$parameters$pro = prob.temp
    #   dens$parameters$mean = mean.temp
    #   dens$parameters$variance$sigma = sigma.temp
    #   dens$modelName = name.temp
    #   # cat("Null")
    # }
    
    # Save parameters for plot.
    param.save$pro[[i]] = dens$parameters$pro
    param.save$mean[[i]] = dens$parameters$mean
    param.save$sigmasq[[i]] = dens$parameters$variance$sigma
    pred.fit[i,] = dens$density
    list_names[i] <- dens$modelName
    
    # density.temp = dens$density
    # prob.temp = dens$parameters$pro
    # mean.temp = dens$parameters$mean
    # sigma.temp = dens$parameters$variance$sigma
    # name.temp = dens$modelName
  }
  
  # Discrepancy function between the predictive surrogate densities 
  # and the posterior predictive density
  discr = matrix(0,kmax,sample_size)
  for(i in 1:kmax){
    discr[i,] = log(pred.fit[i,])-log(pred.post)
  }
  
  # temp = matrix(0,sample_size,kmax)
  # for(i in 1:sample_size){
  #   for(k in 1:kmax){
  #     temp[i,k] = mean(log(pred.fit[k,]) - log(expec.pred.post[i,]))
  #   }  
  # }
  # 
  # boxplot(temp)
  # discr = t(temp)
  # temp1 = apply(discr,1, function(x) quantile(x,0.975))
  # temp2 = apply(discr,1, function(x) quantile(x,0.025))
  
  # temp1 = apply(discr,1, function(x) max(x[abs(x - mean(x)) < sd(x)]))
  # temp2 = apply(discr,1, function(x) min(x[abs(x - mean(x)) < sd(x)]))
  
  temp1 = apply(discr,1, function(x) mean(x) + sd(x))
  temp2 = apply(discr,1, function(x) mean(x) - sd(x))
  
  list_output <- list()
  data.possum = data.frame(upper = temp1,
                           lower = temp2,
                           avg.possum = apply(discr,1,mean),
                           method = "DPM.mult",
                           model_names = list_names,
                           num.fact = 1:kmax)
  
  
  list_output[[1]] <- data.possum
  list_output[[2]] <- param.save
  list_output[[3]] <- unlist(lapply(dp.mult$labelsChain, function(x) length(unique(x))))
  list_output[[5]] <- pred.post
  list_output[[6]] <- y.pred
  list_output[[7]] <- discr
  list_output[[9]] <- dp.mult$labelsChain
  list_output[[10]] <- dp.mult$clusterLabels
  
  # possum.samp = 1000
  cat("Posterior Sample\n")
  # M = R = sample_size
  if(possum.samp != 0){
    H = possum.samp
    # y.pred.quant = matrix(0,N,H)
    y.pred.quant = array(0,dim = c(H,p.dim,sample_size))
    weigths.temp <- list()
    for(i in 1:sample_size){
      
      G <- length(dp.mult$weightsChain[[i]])
      alpha.temp = dp.mult$alphaChain[i]
      weigths.temp = c(dp.mult$weightsChain[[i]]*(n.obs.dpm.mult/(alpha.temp + n.obs.dpm.mult)),alpha.temp/(alpha.temp + n.obs.dpm.mult))
      # weigths.temp = dp.mult$weightsChain[[i]]
      
      # comp = sample(c(1:(G)),H,prob = weigths.temp, replace = TRUE)
      comp = sample(c(1:(G+1)),H,prob = weigths.temp, replace = TRUE)
      
      for(j in 1:H){
        if(comp[j] == (G+1)){
          Lamb.temp = LaplacesDemon::rinvwishart(nu.prior, solve(T0))
          y.pred.quant[j,,i] = mvrnorm(n = 1, mu = mu.prior,Sigma = (1/kappa0)*Lamb.temp)
        }else{
          mu_vec = dp.mult$clusterParametersChain[[i]][1]$mu[,,comp[j]]
          sigma_mat = dp.mult$clusterParametersChain[[i]][2]$sig[,,comp[j]]
          y.pred.quant[j,,i] = mvrnorm(n = 1, mu = mu_vec, Sigma = sigma_mat)
        }
      }
    }
    list_output[[4]] <- y.pred.quant
  }
  
  if(f.pred == TRUE){
    cat("Posterior Expected density\n")
    if(pred.mod == "simulation"){
      y1.dist <- seq(-1, 12.5, length.out = 500)
      y2.dist <- seq(-1, 12.5, length.out = 500)
    }else if(pred.mod == "old_faith"){
      y1.dist <- seq(-0.5, 6, length.out = 500)
      y2.dist <- seq(35, 100, length.out = 500)
    }
    
    dens.post.1 = matrix(0,1000,500)
    dens.post.2 = matrix(0,1000,500)
    
    for(j in 1:1000){
      comp = length(dp.mult$weightsChain[[j]])
      density.temp.1 = numeric(comp)
      density.temp.2 = numeric(comp)
      prob = dp.mult$weightsChain[[j]]
      for(i in 1:500){
        for(k in 1:comp){
          mu_vec = dp.mult$clusterParametersChain[[j]][1]$mu[,,k]
          sigma_mat = dp.mult$clusterParametersChain[[j]][2]$sig[,,k]
          
          density.temp.1[k] = prob[k]*dnorm(y1.dist[i],mean = mu_vec[1],sd = sqrt(sigma_mat[1,1]))
          density.temp.2[k] = prob[k]*dnorm(y2.dist[i],mean = mu_vec[2],sd = sqrt(sigma_mat[2,2]))
          
        }
        dens.post.1[j,i] = sum(density.temp.1)
        dens.post.2[j,i] = sum(density.temp.2)
      }
    }
    expec.dens.post.1 = apply(dens.post.1,2,mean)
    expec.dens.post.2 = apply(dens.post.2,2,mean)
  }
  list_output[[8]] <- cbind(expec.dens.post.1,expec.dens.post.2)
  
  return(list_output)
  
}

#----------------------------------------##----------------------------------------#
# Density summarization - model with dimension d >= 2
#----------------------------------------##----------------------------------------#
# Model-based clustering based on sparse finite Gaussian mixtures - 
# Gertraud Malsiner-Walli - Sylvia Frühwirth-Schnatter - Bettina Grun 2016
# Sparse finite Gaussian mixture model
####################################################################################

# # Multivariate data with cluster information
# y.data = y.data.3
# 
# # columns where the data is located - vector with column numbers
# # col.data
# col.data = 1:2
# 
# # columns where the cluster information is located - vector with column numbers
# # clust.info = c(3)
# 
# clust.info = 0
# y.data = y.data.2
# col.data = 1:2
# kmax = 10
# p.dim = 2
# remove p.dim
# 
# y.data = y.data.2
dcpossum.SFM.mult = function(y.data, kmax, col.data, clust.info = 0, 
                             quant.sample = 1000, possum.samp = 0, 
                             post.size = 1000, f.pred = FALSE, pred.mod = "simulation"){
  
  p.dim = length(col.data)
  
  setwd("/Users/hbolfarine/")
  
  source("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/Sp_Mix_source/Analysis_SpMix_mod_clust.R")
  
  post.list = spmix.mult(y.data, col.data, kmax, clust.info, post.size = post.size)
  
  # Mixture weights
  Eta = post.list$Eta
  
  # Mixture locations
  Mu = post.list$Mu
  
  # Mixture Sigma
  Sigma = post.list$Sigma
  
  # Size of posterior samples
  M = dim(Eta)[1]
  
  # Sample size of the posterior predictive sample
  # I this case will be the same size as the posterior sample
  sample_size <- M
  
  # Posterior predictive sample
  y.pred = matrix(0,sample_size,p.dim)
  for(i in 1:sample_size){
    comp <- sample(1:kmax,1,prob = Eta[i,])
    y.pred[i,] <- mvrnorm(n = 1, mu = Mu[i,,comp], Sigma = Sigma[i,,,comp])
  }
  
  # pairs(y.pred)
  
  # sample_size = 2
  # density.temp <- numeric(kmax)
  # expec.pred.post = matrix(0,sample_size,sample_size)
  # for(j in 1:sample_size){
  #   for(i in 1:sample_size){
  #     for(k in 1:kmax){
  #       density.temp[k] <- Eta[j,k]*dmvnorm(y.pred[i,],mean = Mu[j,,k],sigma = Sigma[j,,,k])
  #     }
  #     expec.pred.post[j,i] = sum(density.temp)
  #   }
  #   # cat(round(j/R,3),"\n")
  #   cat(j)
  # }
  
  setwd("/Users/hbolfarine/")
  sourceCpp("Dropbox/DSS_MIX/DSS_MIX_06_2024/paper_18_07_2024/Rcpp_DSS_Mix/pred_SFM_mult_12_23_2024.cpp")
  cat("cpp-start")
  
  Sigma_mod <- array(Sigma, dim = c(post.size, p.dim*p.dim,kmax))
  
  # Sigma_mod[1,,1]
  # Sigma[1,,,1]
  # sample_size = 1000
  
  # sample_size = 2
  
  expec.pred.post = expect_prob_SFM(sample_size, kmax, p.dim, y.pred, Eta, Mu, Sigma_mod)
  
  pred.post = apply(expec.pred.post,2,mean)
  cat("-end_Pred-")
  #----------------------------------------#
  # Posterior Summaries 
  sum.fit <- numeric(sample_size)
  y.pred.data <- data.frame(y.pred)
  pred.fit = matrix(0,kmax,sample_size)
  list_names <- numeric(kmax)
  param.save = list(prob = list(), mean = list(), sigmasq = list())
  for(i in 1:kmax){
    dens <- densityMclust(y.pred.data, G = i, plot = F, modelNames = "VVV")
    # dens <- densityMclust(y.pred.data, G = i, plot = F)
    
    mod = 1
    while(is.null(dens)){
      dens <- Mclust(y.pred.data, G = i,plot = F, verbose = F, modelNames = model_names[mod])
      mod = mod + 1
    }
    
    # if(length(dens$density) == 0){
    #   dens$density = density.temp
    #   dens$parameters$pro = prob.temp
    #   dens$parameters$mean = mean.temp
    #   dens$parameters$variance$sigma = sigma.temp
    #   dens$modelName = name.temp
    #   # cat("Null")
    # }
    
    # Save parameters for plot.
    param.save$pro[[i]] = dens$parameters$pro
    param.save$mean[[i]] = dens$parameters$mean
    param.save$sigmasq[[i]] = dens$parameters$variance$sigma
    pred.fit[i,] = dens$density
    list_names[i] <- dens$modelName
    
    # density.temp = dens$density
    # prob.temp = dens$parameters$pro
    # mean.temp = dens$parameters$mean
    # sigma.temp = dens$parameters$variance$sigma
    # name.temp = dens$modelName
  }
  
  # Discrenpancy function
  discr = matrix(0,kmax,sample_size)
  for(i in 1:kmax){
    discr[i,] = log(pred.fit[i,])-log(pred.post)
  }
  # 
  
  # temp = matrix(0,M,kmax)
  # for(i in 1:sample_size){
  #   for(k in 1:kmax){
  #     temp[i,k] = mean(log(pred.fit[k,]) - log(expec.pred.post[i,]))
  #   }  
  # }
  
  # boxplot(temp)
  # discr = t(temp)
  # temp1 = apply(discr,1, function(x) quantile(x,0.975))
  # temp2 = apply(discr,1, function(x) quantile(x,0.025))
  
  # temp1 = apply(discr,1, function(x) max(x[abs(x - mean(x)) < sd(x)]))
  # temp2 = apply(discr,1, function(x) min(x[abs(x - mean(x)) < sd(x)]))
  
  temp1 = apply(discr,1, function(x) mean(x) + sd(x))
  temp2 = apply(discr,1, function(x) mean(x) - sd(x))
  
  list_output <- list()
  data.possum = data.frame(upper = temp1,
                           lower = temp2,
                           avg.possum = apply(discr,1,mean),
                           method = "SFM.mult",
                           model_names = list_names,
                           num.fact = 1:kmax)
  
  
  list_output[[1]] <- data.possum
  list_output[[2]] <- param.save
  list_output[[3]] <- post.list$n.comp
  list_output[[5]] <- pred.post
  list_output[[6]] <- y.pred.data
  list_output[[7]] <- discr
  list_output[[9]] <- post.list$clust.SFM
  
  
  # This function relates to the posterior and posterior size M
  if(quant.sample != 0){
    H = quant.sample
    y.pred.quant = array(0,dim = c(H,p.dim,M))
    for(i in 1:M){
      comp = sample(1:kmax,H,prob = Eta[i,],replace = T)
      for(j in 1:H){
        y.pred.quant[j,,i] = mvrnorm(n = 1, mu = Mu[i,,comp[j]], Sigma = Sigma[i,,,comp[j]])
      }
    }
    list_output[[4]] <- y.pred.quant
  }
  cat("-Sample-")
  
  if(f.pred == TRUE){
    
    if(pred.mod == "simulation"){
      y1.dist <- seq(-1, 12.5, length.out = 500)
      y2.dist <- seq(-1, 12.5, length.out = 500)
    }else if(pred.mod == "old_faith"){
      y1.dist <- seq(-0.5, 6, length.out = 500)
      y2.dist <- seq(35, 100, length.out = 500)
    }
    
    dens.post.1 = matrix(0,1000,500)
    dens.post.2 = matrix(0,1000,500)
    
    for(j in 1:1000){
      density.temp.1 = numeric(kmax)
      density.temp.2 = numeric(kmax)
      
      for(i in 1:500){
        for(k in 1:kmax){
          mu_vec = Mu[j,,k]
          sigma_mat = Sigma[j,,,k]
          prob = Eta[j,k]
          
          density.temp.1[k] = prob*dnorm(y1.dist[i],mean = mu_vec[1],sd = sqrt(sigma_mat[1,1]))
          density.temp.2[k] = prob*dnorm(y2.dist[i],mean = mu_vec[2],sd = sqrt(sigma_mat[2,2]))
        }
        dens.post.1[j,i] = sum(density.temp.1)
        dens.post.2[j,i] = sum(density.temp.2)
      }
      
      expec.dens.post.1 = apply(dens.post.1,2,mean)
      expec.dens.post.2 = apply(dens.post.2,2,mean)
    }
    list_output[[8]] <- cbind(expec.dens.post.1,expec.dens.post.2)
  }
  # plot(y2.dist,expec.dens.post.2)
  return(list_output)
  
}

#---------------------------------------------------------------------------#
# y.data = enzyme
# k0 = 0.1
dcpossum.DPM.dir = function(y.data, kmax, quant.sample = 0, k0 = 1, pred.f = TRUE){
  
  # Number of observations 
  n.obs.dpm = length(y.data)
  
  dp.gauss <- DirichletProcessGaussian(y.data, g0Priors = c(median(y.data), k0, 2, 1), alphaPriors = c(2, 4))
  
  # The number controls the number of iterations - 1000
  its = 6000
  dp.inf <- Fit(dp.gauss, its, progressBar = TRUE)
  dp.inf = Burn(dp.inf, 1000)
  # plot(dp.inf)
  
  # Generate posterior predictive sample
  
  R = its-1000
  sample_size = 2000
  y.pred = numeric(sample_size)
  for(i in 1:sample_size){
    G <- length(dp.inf$weightsChain[[i]])
    alpha.temp = dp.inf$alphaChain[i]
    weigths.temp = c(dp.inf$weightsChain[[i]]*(n.obs.dpm/(alpha.temp + n.obs.dpm)),alpha.temp/(alpha.temp + n.obs.dpm))
    
    comp = sample(c(1:(G+1)),1,prob = weigths.temp)
    
    if(comp == (G+1)){
      # Remember to change the base measure
      inv = rinvgamma(1, 2, scale = 1)
      y.pred[i] = rnorm(1,median(y.data),sd = sqrt(inv/k0))
    }else{
      mu = unlist(dp.inf$clusterParametersChain[[i]][1])[comp]
      sigma = unlist(dp.inf$clusterParametersChain[[i]][2])[comp]
      y.pred[i] = rnorm(1,mu,sigma)
    }
  }
  
  expec.pred.post = matrix(0,sample_size,sample_size)
  for(j in 1:sample_size){
    for(i in 1:length(y.pred)){
      mu = unlist(dp.inf$clusterParametersChain[[j]][1])
      sigma = unlist(dp.inf$clusterParametersChain[[j]][2])
      expec.pred.post[j,i] = sum(dp.inf$weightsChain[[j]]*dnorm(y.pred[i],mu,sigma))
    }
    cat(round(j/sample_size,2),"\n")
  }
  
  # Predictive posterior estimate
  pred.post = apply(expec.pred.post,2,mean)
  
  # Generate posterior summaries from 1 to kmax components
  pred.fit = matrix(0,kmax,sample_size)
  param.save = list(prob = list(), mean = list(), sigmasq = list())
  for(i in 1:kmax){
    dens <- densityMclust(y.pred, G = i, plot =F, verbose = FALSE, modelNames = "V")
    # dens <- densityMclust(y.pred, G = i, plot =T, verbose = FALSE)
    dens1 = dens$density
    
    # Save parameters for plot.
    param.save$pro[[i]] = dens$parameters$pro
    param.save$mean[[i]] = dens$parameters$mean
    param.save$sigmasq[[i]] = dens$parameters$variance$sigmasq
    pred.fit[i,] = dens1
  }
  
  # Discrepancy function between the predictive surrogate densities 
  # and the posterior predictive desinty
  discr = matrix(0,kmax,sample_size)
  for(i in 1:kmax){
    discr[i,] = log(pred.fit[i,])-log(pred.post)
  }
  
  # temp = matrix(0,sample_size,kmax)
  # for(i in 1:sample_size){
  #   for(k in 1:kmax){
  #     temp[i,k] = mean(log(pred.fit[k,]) - log(expec.pred.post[i,])) #- mean(log(pred.post))
  #   }  
  # }
  
  
  # log.pred = log(expec.pred.post)
  # discr = matrix(0,kmax,sample_size)
  # log.pred.fit = log(pred.fit)
  # for(i in 1:kmax){
  #   temp.mat = t(log.pred.fit[i,] - t(log.pred))
  #   discr[i,] = apply(temp.mat,1,mean)
  # }
  # 
  # discr = t(temp)
  # boxplot(t(discr))
  # avg.possum = apply(discr,1,mean)
  # plot(avg.possum)
  # temp1 = apply(discr,1, function(x) quantile(x,0.975))
  # temp2 = apply(discr,1, function(x) quantile(x,0.025))
  
  # temp1 = apply(discr,1, function(x) max(x[abs(x - mean(x)) < sd(x)]))
  # temp2 = apply(discr,1, function(x) min(x[abs(x - mean(x)) < sd(x)]))
  
  temp1 = apply(discr,1, function(x) mean(x) + sd(x))
  temp2 = apply(discr,1, function(x) mean(x) - sd(x))
  
  # Saving outputs from the density estimation
  #-------------------------------------------------#
  list_output <- list()
  data.possum = data.frame(upper = temp1,
                           lower = temp2,
                           avg.possum = apply(discr,1,mean),
                           method = "DPM",
                           num.fact = 1:kmax)
  
  
  list_output[[1]] <- data.possum
  list_output[[2]] <- param.save
  list_output[[3]] <- unlist(lapply(dp.inf$labelsChain, function(x) length(unique(x))))
  list_output[[5]] <- pred.post
  list_output[[6]] <- y.pred
  list_output[[7]] <- discr
  
  #-------------------------------------------------# 
  
  if(quant.sample != 0){
    H = quant.sample
    y.pred.quant = matrix(0,sample_size,H)
    weight_list <- list()
    for(i in 1:sample_size){
      
      G <- length(dp.inf$weightsChain[[i]])
      alpha.temp = dp.inf$alphaChain[i]
      weigths.temp = c(dp.inf$weightsChain[[i]]*(n.obs.dpm/(alpha.temp + n.obs.dpm)),alpha.temp/(alpha.temp + n.obs.dpm))
      
      comp = sample(c(1:(G+1)),H,prob = weigths.temp, replace = TRUE)
      
      for(j in 1:H){
        if(comp[j] == (G+1)){
          inv = rinvgamma(1, 2, scale = 1)
          y.pred.quant[i,j] = rnorm(1,median(y.data),sd = sqrt(inv/k0))
        }else{
          mu = unlist(dp.inf$clusterParametersChain[[i]][1])[comp[j]]
          sigma = unlist(dp.inf$clusterParametersChain[[i]][2])[comp[j]]
          y.pred.quant[i,j]  = rnorm(1,mu,sigma)
        }
      }
    }
    list_output[[4]] <- y.pred.quant
  }
  
  if(pred.f == TRUE){
    min.y<-min(y.data)-0.5*sqrt(var(y.data))
    max.y<-max(y.data)+0.5*sqrt(var(y.data))
    y.seq = seq(min.y, max.y, length.out = 500)
    cat("Pred_F")
    expec.pred.post.f = matrix(0,sample_size,500)
    for(j in 1:sample_size){
      # cat(j)
      for(i in 1:500){
        # mat_mu <- param_list[[j]]$theta_mat[,1]
        # mat.sigma <- param_list[[j]]$theta_mat[,2]
        # expec.pred.post.f[j,i]= sum(param_list[[j]]$vector_weights*dnorm(y.seq[i],mean = mat_mu,sd = mat.sigma))
        mu = unlist(dp.inf$clusterParametersChain[[j]][1])
        sigma = unlist(dp.inf$clusterParametersChain[[j]][2])
        expec.pred.post.f[j,i] = sum(dp.inf$weightsChain[[j]]*dnorm(y.seq[i],mu,sigma))
      }
    }
    # plot(y.seq,pred.func.f)
    # Predictive posterior estimate
    pred.func.f = apply(expec.pred.post.f,2,mean)
    list_output[[8]] <- pred.func.f
  }
  
  # expected predictive density - Hellinger
  
  expec.pred.post.data = matrix(0,length(dp.inf$clusterParametersChai),length(y.data))
  for(j in 1:length(dp.inf$clusterParametersChai)){
    # cat(j)
    for(i in 1:length(y.data)){
      # mat_mu <- param_list[[j]]$theta_mat[,1]
      # mat.sigma <- param_list[[j]]$theta_mat[,2]
      # expec.pred.post.f[j,i]= sum(param_list[[j]]$vector_weights*dnorm(y.seq[i],mean = mat_mu,sd = mat.sigma))
      mu = unlist(dp.inf$clusterParametersChain[[j]][1])
      sigma = unlist(dp.inf$clusterParametersChain[[j]][2])
      expec.pred.post.data[j,i] = sum(dp.inf$weightsChain[[j]]*dnorm(y.data[i],mu,sigma))
    }
  }
  # plot(y.seq,pred.func.f)
  # Predictive posterior estimate
  pred.func.expec = apply(expec.pred.post.data,2,mean)
  list_output[[9]] <- pred.func.expec
  
  
  return(list_output)
  # plot(y.seq,pred.func.f)
}


