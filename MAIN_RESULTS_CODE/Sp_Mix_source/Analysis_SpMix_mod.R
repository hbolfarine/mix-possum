###################################################################################
#          Estimation of a sparse finite Gaussian mixture model - main file
####################################################################################

## read sources
library(e1071)
library(mclust)
library(MASS)
library(bayesm)
library(MCMCpack)
library(mvtnorm)
library(Runuran)
library(flexclust)


# MCMC estimation
setwd("/Users/hb23255/")
source("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/Sp_Mix_source/Estimation_SpMix.R")
source("Dropbox/DSS_MIX/DSS_MIX_06_2024/source/Sp_Mix_source/Identification_SpMix.R")

  
spmix.mult = function(y.data, col.data, kmax, clust.info = 0, post.size = 2000){
  
  # Matrix with data
  y <- y.data[,col.data]
  
  # Dimension of the data
  r <- length(y[1, ])
  
  # Number of observations
  N <- length(y[, 1])
  
  # Create variable names
  names_y <- c()
  for (j in 1:r) {
    names_y <- c(names_y, paste("y", j, sep = ""))
  }
  colnames(y) <- names_y
  
  # Cluster tags
  # if(!clust.info == 0){
  #   z <- y.data[,3]
  # }
  
  # Create a numerical matrix 
  # this will be the input for the method 
  y <- as.matrix(y)
  
  #---------- B) Specification of the simulation and prior parameters -----------------------------------------------
  
  ## number of mixture components 
  K <- kmax
  ## number of iterations, M without burnin
  M <- post.size
  burnin <- 200
  
  # Number of iterations with burnin
  Mmax <- M + burnin
  
  ## Dirichlet parameter for the mixture weights
  e0 <- 0.01
  ## variance of the normal proposal for the MH step for estimating e0
  # I don't know what this is?
  c_proposal <- 0.8
  
  # Vector containing the minimum and maximum for each column
  # This value will be used on the prior for the covariance matrix
  # and on the mean vector 
  R <- apply(y, 2, function(x) diff(range(x)))
  
  ## prior on Sigma_k
  c0 <- 2.5 + (r - 1)/2
  C0 <- 0.75 * cov(y)
  g0 <- 0.5 + (r - 1)/2
  G0 <- 100 * g0/c0 * diag((1/R^2))
  
  ## prior on mu
  b0 <- apply(y, 2, median)
  B_0 <- rep(1, r)  #initial values for lambda are 1
  cat(B_0)
  B0 <- diag((R^2) * B_0)
  nu <- 0.5
  
  
  ## initial values for parameters to be estimated:
  eta_0 <- rep(1/K, K)
  sigma_0 <- array(0, dim = c(r, r, K))
  for (k in 1:K) {
    sigma_0[, , k] <- C0
  }
  C0_0 <- C0
  
  ## initial classification
  groups <- K
  cl_y <- kmeans(y, centers = groups, nstart = 30)
  S_0 <- cl_y$cluster
  mu_0 <- cbind(t(cl_y$centers))
  
  
  ## generate matrices for saving the results:
  
  # posteriori for eta 
  Eta <- matrix(0, M, K)
  
  # posteriori for mu
  Mu <- array(0, dim = c(M, r, K))
  
  # posteriori for B
  B <- matrix(0, M, r)
  
  # What are these? - FS?
  # Eta_Matrix_FS <- matrix(0, r, K)
  # Sigma_Matrix_FS <- array(0, dim = c(r, r, K))
  # Mu_Matrix_FS <- array(0, dim = c(r, K))
  
  
  #---------- C) Gibbs sampling from the posterior -----------------------------------------------
  # cat("begin")
  ################ call MCMC procedure
  estGibbs <- MultVar_NormMixt_Gibbs_IndPriorNormalgamma(y, S_0, mu_0, sigma_0, eta_0, e0, c0, C0_0, 
                                                         g0, G0, b0, B0, nu, B_0, M, burnin, c_proposal, 
                                                         priorOnE0 = FALSE, lambda = FALSE, R)
  # cat("end")
  Mu <- estGibbs$Mu
  Sigma <- estGibbs$Sigma
  Eta <- estGibbs$Eta
  # S_alt_matrix <- estGibbs$S_alt_matrix
  # B <- estGibbs$B
  # e0_vector <- estGibbs$e0_vector
  # acc_rate <- estGibbs$acc_rate
  Nk_matrix_alt <- estGibbs$Nk_matrix_alt
  # nonnormpost_mode_list <- estGibbs$nonnormpost_mode_list
  
  #------------------------------------------#
  K0_vector <- rowSums(Nk_matrix_alt != 0)  #vector with number of non-empty groups of each iteration
  # p_K0 <- tabulate(K0_vector, K)
  # p_K0
  # par(mfrow = c(1, 1))
  # barplot(p_K0, names = 1:K, xlab = "number of non-empty groups K0", col = "green", ylab = "freq")
  # K0 <- which.max(p_K0)
  # K0  #mode K0 is the estimator for K_true
  # M0 <- sum(K0_vector == K0)
  # M0  #M0 draws have exactly K0 non-empty groups
  
  # 
  list.SFM = list(Eta = Eta, Mu = Mu, Sigma = Sigma, n.comp = K0_vector)
  
  return(list.SFM)
}
#------------------------------------------------#

# M = 1000
# G = 10
# p.dim = dim(y)[2]
# # Posterior samples from mean bayesm
# # Posterior predictive sample
# y.pred = matrix(0,M,p.dim)
# for(i in 1:M){
#   comp <- sample(1:G,1,prob = Eta[i,])
#   y.pred[i,] <- mvrnorm(n = 1, mu = Mu[i,,comp], Sigma = Sigma[i,,,comp])
# }
# 
# pairs(y.pred)
# 
# density.temp <- numeric(G)
# expec.pred.post = matrix(0,M,M)
# for(j in 1:M){
#   for(i in 1:M){
#     for(k in 1:G){
#       density.temp[k] <- Eta[j,k]*dmvnorm(y.pred[i,],mean = Mu[j,,k],sigma = Sigma[j,,,k])
#     }
#     expec.pred.post[j,i] = sum(density.temp)
#   }
#   # cat(round(j/R,3),"\n")
#   cat(j)
# }

# 
# pred.post = apply(expec.pred.post,2,mean)
# 
# #----------------------------------------#
# kmax = 10
# # Posterior Summaries 
# sum.fit <- numeric(M)
# y.pred.data <- data.frame(y.pred)
# pred.fit = matrix(0,kmax,M)
# list_names <- numeric(kmax)
# param.save = list(prob = list(), mean = list(), sigmasq = list())
# for(i in 1:kmax){
#   dens <- densityMclust(y.pred.data, G = i,plot = F)
#   dens1 = dens$density
#   
#   # Save parameters for plot.
#   param.save$pro[[i]] = dens$parameters$pro
#   param.save$mean[[i]] = dens$parameters$mean
#   param.save$sigmasq[[i]] = dens$parameters$variance$sigma
#   pred.fit[i,] = dens1
#   list_names[i] <- dens$modelName
# }
# 
# # Discrenpancy function
# discr = matrix(0,kmax,M)
# for(i in 1:kmax){
#   discr[i,] = log(pred.fit[i,])-log(pred.post)
# }
# 
# temp1 = apply(discr,1, function(x) quantile(x,0.975))
# temp2 = apply(discr,1, function(x) quantile(x,0.025))
# 
# # temp1 = apply(discr,1, function(x) mean(x) + sd(x))
# # temp2 = apply(discr,1, function(x) mean(x) - sd(x))
# 
# temp.22 = discr
# sel.K = max(which(apply(temp.22,1,mean) < temp2[G]))+1
# # sel.K
# diff.sd = temp1-temp2
# sel.K = max(which(apply(temp.22,1,mean) < temp2[which(diff.sd == min(diff.sd))]))+1
# sel.K
# 
# plot(apply(temp.22,1,mean), type = "l", ylim = c(min(temp2)+0.1*min(temp2), max(temp1)+0.1*max(temp1)), main = "Discrepancy plot", ylab = "Discrepancy", xlab = "Number of components")
# for(i in 1:G){
#   segments(i,temp1[i],i,temp2[i])
# }
# points(apply(temp.22,1,mean))
# abline(h = q1, col = "blue", lty = 2)
# abline(h = 0, col = "red")
# abline(h = temp2[G], col = "green", lty = 2)
# 
