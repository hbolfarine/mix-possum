library(dirichletprocess)
library(invgamma)
library(LaplacesDemon)

dpm.dir.mult = function(y.data, kappa0 = 1, nu.prior = 2){
  
  # old, kappa = 10, nu = 5
  # old, kappa = 1, nu = 2
  # sim, kappa = 1, nu = 2
  # thyroid, kappa = 1, nu = 5
  
  # Number of observations in the data
  n.obs.dpm.mult = dim(y.data)[1]
  
  # Prior mean 
  mu.prior = apply(y.data, 2, mean)
  
  # Covariance prior
  # T0 = diag(ncol(y.data.2))
  T0 = cov(y.data)
  
  # kappa0 = 1
  # # kappa0 = 2
  # nu.prior = 5
  # # nu.prior = 2

  g0Priors <- list(mu0 = mu.prior,
                   Lambda = T0,
                   kappa0 = kappa0,
                   nu = nu.prior) 
  
  dp.mult <- DirichletProcessMvnormal(as.matrix(y.data.2), g0Priors = g0Priors,
                                      alphaPriors = c(2,4))
  
  its = 2000
  dp.mult <- Fit(dp.mult, its, progressBar = TRUE)
  dp.mult = Burn(dp.mult, 1000)
  plot(dp.mult)
  
}