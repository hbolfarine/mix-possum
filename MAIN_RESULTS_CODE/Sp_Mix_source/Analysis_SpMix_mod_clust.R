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
setwd("/Users/hbolfarine/")
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
  S_alt_matrix <- estGibbs$S_alt_matrix
  # B <- estGibbs$B
  # e0_vector <- estGibbs$e0_vector
  # acc_rate <- estGibbs$acc_rate
  Nk_matrix_alt <- estGibbs$Nk_matrix_alt
  nonnormpost_mode_list <- estGibbs$nonnormpost_mode_list
  
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
  # list.SFM = list(Eta = Eta, Mu = Mu, Sigma = Sigma, n.comp = K0_vector)
  
  #------------------------------------------------#
  
  
  ##### number of nonempty components (nne_gr), after burnin:
  # nne_gr <- apply(Nk_matrix_alt, 1, function(x) sum(x != 0))
  # table(nne_gr)
  # plot(1:length(nne_gr), nne_gr, type = "l", ylim = c(1, max(nne_gr)), 
  #      xlab = paste("iteration"), main = "number of non-empty components")
  
  
  
  #---------- E) Identification of the mixture model -----------------------------------------------
  
  ##### estimating the number of nonempty components
  K0_vector <- rowSums(Nk_matrix_alt != 0)  #vector with number of non-empty groups of each iteration
  p_K0 <- tabulate(K0_vector, K)
  p_K0
  # par(mfrow = c(1, 1))
  # barplot(p_K0, names = 1:K, xlab = "number of non-empty groups K0", col = "green", ylab = "freq")
  K0 <- which.max(p_K0)
  K0  #mode K0 is the estimator for K_true
  M0 <- sum(K0_vector == K0)
  M0  #M0 draws have exactly K0 non-empty groups
  
  
  ##### selecting those draws where the number of non-empty groups was exactly K0:
  Nk_matrix_K0 <- (Nk_matrix_alt[K0_vector == K0, ] != 0)
  M0_Nk_1 <- sum(rowSums((Nk_matrix_alt[K0_vector == K0, ]) == 1) == 1)
  M0_Nk_1  ##how many draws have components with one observation only?
  
  
  ##### extracting those draws which are sampled from exactly K0 non-empty groups 
  #### i)Mu:
  Mu_inter <- array(0, dim = c(M0, r, K))
  Mu_inter <- Mu[K0_vector == K0, , ]  #matrix with draws only from the K0 interesting groups
  Mu_K0 <- array(0, dim = c(M0, r, K0))  #matrix where empty groups are cut off
  for (j in 1:r) {
    Mu_K0[, j, ] <- matrix(t(Mu_inter[, j, ])[t(Nk_matrix_K0)], ncol = K0, byrow = T)
  }
  
  
  #### iii)Eta:
  Eta_inter <- matrix(0, M0, K)  #matrix with draws only from the K0 interesting groups
  Eta_inter <- Eta[K0_vector == K0, ]
  Eta_K0 <- matrix(t(Eta_inter)[t(Nk_matrix_K0)], byrow = T, ncol = K0)
  
  
  #### iv)S_matrix:
  S_matrix_inter <- matrix(0, M, N)
  for (m in 1:M) {
    if (K0_vector[m] == K0) {
      perm_S <- rep(0, K)
      perm_S[Nk_matrix_alt[m, ] != 0] <- 1:K0
      S_matrix_inter[m, ] <- perm_S[S_alt_matrix[m, ]]
    }
  }
  S_matrix_K0 <- S_matrix_inter[K0_vector == K0, ]
  
  
  ##### clustering and relabeling of the MCMC draws corresponding to the selected model
  map_mu <- nonnormpost_mode_list[[K0]]$bk  #initial values for the group centers
  map_cov <- nonnormpost_mode_list[[K0]]$Bk  #initial values for the group covariances
  clust_FS_K0 <- MultVar_clust_FS_FAST_2(Mu_K0, Eta_K0, S_matrix_K0, map_mu, map_cov, maha = FALSE)
  Mu_only_perm <- clust_FS_K0$Mu_only_perm
  Eta_only_perm <- clust_FS_K0$Eta_only_perm
  S_matrix_only_perm <- clust_FS_K0$S_matrix_only_perm
  non_perm_rate <- clust_FS_K0$non_perm_rate
  # non_perm_rate         #non-permuation rate
  Mu_Matrix_FS <- colMeans(Mu_only_perm, dims = 1)  #estimated component means
  Eta_Matrix_FS <- colMeans(Eta_only_perm, dims = 1)  #estimated mixture weights
  
  
  
  #---------- F) Classification of the observations -----------------------------------------------
  
  ##### frequency of the assignment to the groups for each observation:
  Ass <- matrix(0, N, K0)
  for (n in 1:N) {
    Ass[n, ] <- tabulate(S_matrix_only_perm[, n], K0)
  }
  ass <- c()
  ass <- apply(Ass, 1, which.max)  #estimed classification
  
  list.SFM = list(Eta = Eta, Mu = Mu, Sigma = Sigma, n.comp = K0_vector, clust.SFM = ass)
  
  return(list.SFM)
}
##### evaluate the estimated classification
# table(z, ass)
# classError(ass, z)$errorRate
# adjustedRandIndex(z, ass)
