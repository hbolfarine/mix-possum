#----------------------------------------#
# Data for simulation
#----------------------------------------#

# n = 1000 generated from the eighth density function in Table 1 of Marron & Wand (1992):
library(densEstBayes)
# x <- rMarronWand(1000,8)

# https://search.r-project.org/CRAN/refmans/nor1mix/html/MarronWand.html

data_sim_func = function(data_obs,n = 1000){
  
  if(data_obs == "wand"){
    sample = numeric(n)
    
    sample = rMarronWand(1000,8)
    
    return(sample)
    
  }else if(data_obs == "example_01"){
    sample = numeric(n)
    
    for(i in 1:n){
      u = runif(1)
      if(u <= 0.5){
        sample[i] = rnorm(1,1,sqrt(1))
      }else if(u > 0.5){
        sample[i] = rnorm(1,8,sqrt(1))
      }
    }
    
    return(sample)
    
  }else if(data_obs == "example_02"){
    sample = numeric(n)
    
    # 0.2 \times N(19,5) + 0.2 \times N(19,1) + 0.25\times N(23,1) + 0.2\times N(29,1) + 0.15\times N(33,2)
    
    for(i in 1:n){
      u = runif(1)
      if(u < 0.2){
        sample[i] = rnorm(1,19,sqrt(5))
      }else if(u >= 0.2 &&  u < 0.4){
        sample[i] = rnorm(1,19,sqrt(1))
      }else if(u >= 0.4 &&  u < 0.65){
        sample[i] = rnorm(1,23,sqrt(1))
      }else if(u >= 0.65 &&  u < 0.85){
        sample[i] = rnorm(1,29,sqrt(1))
      }else if(u >= 0.85){
        sample[i] = rnorm(1,33,sqrt(2))
      }
    }
    
    # for(i in 1:n){
    #   u = runif(1)
    #   if(u < 0.2){
    #     sample[i] = rnorm(1,10,sqrt(10))
    #   }else if(u >= 0.2 &&  u < 0.8){
    #     sample[i] = rnorm(1,19,sqrt(2))
    #   }else if(u >= 0.8){
    #     sample[i] = rnorm(1,23,sqrt(1))
    #   }
    # }
    
    return(sample)
    # 0.2×N (1, 1)+0.6×N (3, 6)+ 0.2×N(10,2)
  }else if(data_obs == "example_03"){
    sample = numeric(n)
    
    for(i in 1:n){
      u = runif(1)
      if(u < 0.2){
        sample[i] = rnorm(1,1,sqrt(1))
      }else if(u >= 0.2 &&  u < 0.8){
        sample[i] = rnorm(1,3,sqrt(6))
      }else if(u >= 0.8){
        sample[i] = rnorm(1,10,sqrt(2))
      }
    }
    
    return(sample)
    
  }else if(data_obs == "examp_laplace"){
    
    sample = numeric(n)
    # means (−5, 5), scales (1.5, 1), and mixing weights (0.4, 0.6). 
    for(i in 1:n){
      u = runif(1)
      if(u <= 0.4){
        sample[i] = LaplacesDemon::rlaplace(1, location= -5, scale = 1.5)
      }else if(u > 0.4){
        sample[i] = LaplacesDemon::rlaplace(1, location= 5, scale = 1)
      }
    }
    
    return(sample)  
    
  }else if(data_obs == "example_mult_02"){
    
    sample = matrix(0,ncol = 2, nrow = n)
    
    for(i in 1:n){
      u = runif(1)
      if(u <= 0.5){
        mean_vector = c(10, 10)
        covariance_matrix <- matrix(c(1, 0, 0, 1), ncol = 2)
      }else{
        mean_vector = c(0, 5)
        covariance_matrix <- matrix(c(1, 0, 0, 1), ncol = 2)
      }
      sample[i,] <- mvrnorm(n = 1, mu = mean_vector, Sigma = covariance_matrix)
    }
    
    return(sample)
    
  }else if(data_obs == "galaxy"){  
    
    data("galaxies", package = "MASS")
    # now fix a typographical error in the data
    # see help("galaxies", package = "MASS")
    galaxies[78] = 26960
    appl.data = galaxies / 1000
    
    return(appl.data)
    
  }else if(data_obs == "acidity"){
    
    data(acidity, package = "mclust")
    appl.data = acidity
    return(appl.data)
    
  }else if(data_obs == "faithful"){
    
    # data(faithful, package = "mvtnorm")
    appl.data <- data.frame(faithful)
    # appl.data <- as.matrix(appl.data)
    return(appl.data)
    
  }else if(data_obs == "thyroid"){
    
    data(thyroid, package = "mclust")
    appl.data <- data.frame(thyroid)
    appl.data <- as.matrix(appl.data[,2:6])
    return(appl.data)
    
  }else if(data_obs == "crabs"){
    
    data(crabs, package = "MASS")
    appl.data <- data.frame(crabs)
    appl.data <- as.matrix(appl.data[,4:8])
    return(appl.data)
    
  }else if(data_obs == "example_mult_03"){
    
    sample = matrix(0,ncol = 2, nrow = n)
    
    for(i in 1:n){
      u = runif(1)
      if(u <= 0.45){
        mean_vector = c(4, 4) 
        covariance_matrix <- matrix(c(1, 0, 0, 1), ncol = 2)
      }else if(0.45<= u & u <= 0.75){
        mean_vector = c(7, 4) 
        covariance_matrix_01 <- matrix(c(2.5, 0, 0, 0.2), ncol = 2)
        Rotate_matrix <- matrix(c(cos(pi/4),sin(pi/4),-sin(pi/4),cos(pi/4)), ncol = 2)
        covariance_matrix <- Rotate_matrix%*%covariance_matrix_01%*%t(Rotate_matrix)
      }else if(u >= 0.75){
        mean_vector = c(6, 2) 
        covariance_matrix <- matrix(c(3, 0, 0, 0.1), ncol = 2)
      }
      sample[i,] <- mvrnorm(n = 1, mu = mean_vector, Sigma = covariance_matrix)
    }
    return(sample)
    
  }else if(data_obs == "example_01_clust"){
    sample = numeric(n)
    clust_assign = numeric(n)
    
    for(i in 1:n){
      u = runif(1)
      if(u <= 0.5){
        sample[i] = rnorm(1,1,sqrt(1))
        clust_assign[i] = 1
      }else if(u > 0.5){
        sample[i] = rnorm(1,8,sqrt(1))
        clust_assign[i] = 2
      }
    }
    return(cbind(sample,clust_assign))
    
  }else if(data_obs == "example_02_clust"){
    sample = numeric(n)
    clust_assign = numeric(n)
    
    # 0.2 \times N(16,5) + 0.2 \times N(19,1) + 0.25\times N(23,1) + 0.2\times N(29,1) + 0.15\times N(33,1)
    
    for(i in 1:n){
      u = runif(1)
      if(u < 0.2){
        sample[i] = rnorm(1,19,sqrt(5))
        clust_assign[i] = 1
      }else if(u >= 0.2 &&  u < 0.4){
        sample[i] = rnorm(1,19,sqrt(1))
        clust_assign[i] = 2
      }else if(u >= 0.4 &&  u < 0.65){
        sample[i] = rnorm(1,23,sqrt(1))
        clust_assign[i] = 3
      }else if(u >= 0.65 &&  u < 0.85){
        sample[i] = rnorm(1,29,sqrt(1))
        clust_assign[i] = 4
      }else if(u >= 0.85){
        sample[i] = rnorm(1,33,sqrt(2))
        clust_assign[i] = 5
      }
    }
    return(cbind(sample,clust_assign))
    
  }else if(data_obs == "example_mult_03_clust"){
    
    sample = matrix(0,ncol = 2, nrow = n)
    clust_assign = numeric(n)
    
    for(i in 1:n){
      u = runif(1)
      if(u <= 0.45){
        mean_vector = c(4, 4) 
        covariance_matrix <- matrix(c(1, 0, 0, 1), ncol = 2)
        clust_assign[i] = 1
      }else if(0.45<= u & u <= 0.75){
        mean_vector = c(7, 4) 
        covariance_matrix_01 <- matrix(c(2.5, 0, 0, 0.2), ncol = 2)
        Rotate_matrix <- matrix(c(cos(pi/4),sin(pi/4),-sin(pi/4),cos(pi/4)), ncol = 2)
        covariance_matrix <- Rotate_matrix%*%covariance_matrix_01%*%t(Rotate_matrix)
        clust_assign[i] = 3
      }else if(u >= 0.75){
        mean_vector = c(6, 2) 
        covariance_matrix <- matrix(c(3, 0, 0, 0.1), ncol = 2)
        clust_assign[i] = 2
      }
      sample[i,] <- mvrnorm(n = 1, mu = mean_vector, Sigma = covariance_matrix)
    }
    return(cbind(sample,clust_assign))
    
  }else if(data_obs == "example_mult_02_clust"){
    
    sample = matrix(0,ncol = 2, nrow = n)
    clust_assign = numeric(n)
    
    for(i in 1:n){
      u = runif(1)
      if(u <= 0.5){
        mean_vector = c(6, 6)
        covariance_matrix <- matrix(c(1, 0, 0, 1), ncol = 2)
        clust_assign[i] = 2
      }else{
        mean_vector = c(5, 5)
        covariance_matrix <- matrix(c(1, 0, 0, 1), ncol = 2)
        clust_assign[i] = 1
      }
      sample[i,] <- mvrnorm(n = 1, mu = mean_vector, Sigma = covariance_matrix)
    }
    return(cbind(sample,clust_assign))
    
  }else if(data_obs == "thyroid_clust"){
    
    data(thyroid, package = "mclust")
    appl.data <- data.frame(thyroid)
    return(appl.data)
    
  }else if(data_obs == "crabs_clust"){
    
    data(crabs, package = "MASS")
    appl.data <- data.frame(crabs)
    return(appl.data)
    
  }

}
