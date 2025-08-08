model_names <- c(
  "VVV",   # ellipsoidal, varying volume, shape, and orientation
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

assign_clusters_mult <- function(data, means) {
  distances <- as.matrix(dist(rbind(means, data)))
  distances <- distances[-seq(nrow(means)), 1:nrow(means)]
  cluster_assignments <- apply(distances, 1, which.min)
  return(cluster_assignments)
}

# data = y.data
# means = fixed_means
assign_clusters <- function(data, means) {
  distances <- abs(outer(data, means[,1], "-"))
  cluster_assignments <- apply(distances, 1, which.min)
  return(cluster_assignments)
}

# list.dc.possum = list_output
# list.dc.possum = BP.acidity
# y.data = y.data.2
# y.data = y.data.1
# make function with predict
dc.possum.clust = function(list.dc.possum, y.data, K.sel, km = FALSE){
  
  y.fit.pred = list.dc.possum[[4]]
  
  # H,p.dim,N
  N = dim(y.fit.pred)[3]
  H = dim(y.fit.pred)[3]
  
  #------------------------------------------------------#
  
  N.samp = nrow(y.data)
  mat.respons = matrix(0,N.samp,K.sel)
  mat.clust = matrix(0,N,N.samp)
  colnames(mat.respons) = c(1:K.sel)
  mat.order = matrix(0,N,K.sel)
  
  #------------------------------------------------------#
  p.dim = dim(y.data)[2]
  if(p.dim == 2){
    data.orig <- data.frame(x = y.data[,1],y = y.data[,2])
    mat.clust.k = matrix(0,N,N.samp)
  }
  
  # predict(clust,y.data)$classification
  
  # comprar com densityMclust
  for(i in 1:N){
    clust <- Mclust(y.fit.pred[,,i], G = K.sel,plot = F, verbose = F, modelNames = "VVV")
    
    mod = 1
    while(is.null(clust)){
      clust <- Mclust(y.fit.pred[,,i], G = K.sel,plot = F, verbose = F, modelNames = model_names[mod])
      mod = mod + 1
    }
    
    order = order(clust$parameters$mean[1,])
    mat.order[i,] = order
    # order = order(xtabs(~predict(clust,y.data)$classification))
    # mat.order[i,] = order
    mean.temp = clust$parameters$mean[,order]
    prob.temp = clust$parameters$pro[order]
    sigma.temp = clust$parameters$variance$sigma[,,order]
    
    # mean.temp = clust$parameters$mean
    # prob.temp = clust$parameters$pro
    # sigma.temp = clust$parameters$variance$sigma
    

    
    for(k in 1:K.sel){
      mat.respons[,k] = prob.temp[k]*dmvnorm(y.data,mean = mean.temp[,k], sigma = sigma.temp[,,k])
    }
    mat.respons = mat.respons/apply(mat.respons, 1, sum)
    
    temp.clust = as.numeric(apply(mat.respons, 1, which.max))
    mat.clust[i,] = temp.clust # important
    # mat.clust[i,] = predict(clust,y.data)$classification
    # cat(i)
    
    if(p.dim == 2){
      kmeans_result <- kmeans(y.fit.pred[,,i], centers = K.sel)
      centers = t(kmeans_result$centers)
      order = order(centers[1,])
      centers = centers[,order]
      fixed_means <- data.frame(
        x = centers[1,],
        y = centers[2,]
      )
      
      mat.clust.k[i,] <- as.numeric(assign_clusters_mult(data.orig, fixed_means))
      cat(i/N)
    }
    cat(i/N)
  }
  
  temp = apply(mat.clust,2,function(x) 1 - max(xtabs(~x)/N))
  data.clust = data.frame(dat = y.data, clust.prob = temp)
  
  if(km == TRUE){
    prob.km = apply(mat.clust.k,2,function(x) 1 - max(xtabs(~x)/N))
    data.clust = cbind(data.clust, clust.prob.km = prob.km)
  }
  
  return(data.clust)
  
}

# list.dc.possum = BP.enzyme
# y.data = enzyme
dc.optim.clust.uni = function(list.dc.possum, y.data, K.sel){
  
  parm = list.dc.possum[[2]]
  mean = parm$mean[[K.sel]]
  prob = parm$pro[[K.sel]]
  sigma.optim = parm$sigmasq[[K.sel]]
  
  # if(sum(is.na(sigma.optim)) != 0){
  #   sigma.optim[1:K.sel] = sigma.temp[which(!is.na(sigma.temp))]
  # }
  
  N.samp = length(y.data)
  mat.respons = matrix(0,N.samp,K.sel)
  
  for(k in 1:K.sel){
    mat.respons[,k] = prob[k]*dnorm(y.data,mean = mean[k], sd = sigma.optim[k])
  }
  
  mat.respons.prob = mat.respons/apply(mat.respons, 1, sum)
  
  optim.clust = apply(mat.respons.prob, 1, which.max)
  
  data.clust = data.frame(dat = y.data, optim.clust = optim.clust)
  
  return(data.clust)
  
}


dc.optim.clust.kmeans = function(list.dc.possum, y.data, K.sel){
  
  y.pred = list.dc.possum[[6]]
  
  data.orig <- data.frame(x = y.data[,1],y = y.data[,2])
  
  N.samp = nrow(y.data)
  mat.respons = matrix(0,N.samp,K.sel)
  
  kmeans_result <- kmeans(y.pred, centers = K.sel)
  centers = t(kmeans_result$centers)
  order = order(centers[1,])
  centers = centers[,order]
  fixed_means <- data.frame(
    x = centers[1,],
    y = centers[2,]
  )
  
  optim.clust.k <- as.numeric(assign_clusters_mult(data.orig, fixed_means))
  
  data.clust = data.frame(dat = y.data, clust.km = optim.clust.k)
  
  return(data.clust)
  
}

# list.dc.possum = DPM.galaxy
# # # # y.data = y.data.2
# y.data = y.data.app
# y.data = acidity
dc.possum.clust.uni = function(list.dc.possum, y.data, K.sel, km = FALSE){
  
  y.fit.pred = list.dc.possum[[4]]
  
  # list.dc.possum[[2]]
  
  # H,p.dim,N
  N = dim(y.fit.pred)[2]
  H = dim(y.fit.pred)[2]
  
  #------------------------------------------------------#
  
  N.samp = length(y.data)
  mat.respons = matrix(0,N.samp,K.sel)
  mat.clust = matrix(0,N,N.samp)
  mat.clust.k = matrix(0,N,N.samp)
  colnames(mat.respons) = c(1:K.sel)
  mat.order = matrix(0,N,K.sel)
  
  #------------------------------------------------------#
  # p.dim = dim(y.data)[2]
  # if(p.dim == 2){
  #   data.orig <- data.frame(x = y.data[,1],y = y.data[,2])
  #   mat.clust.k = matrix(0,N,N.samp)
  # }
  
  par <- vector("list", 3)
  par$pro <- list.dc.possum[[2]]$pro[[K.sel]]
  par$mean <- list.dc.possum[[2]]$mean[[K.sel]]
  # par$sigma <- list.dc.possum[[2]]$sigmasq[[K.sel]]
  par$variance$sigmasq <- list.dc.possum[[2]]$sigmasq[[K.sel]]
  
  # comprar com densityMclust
  for(i in 1:N){
    clust <- Mclust(y.fit.pred[i,], G = K.sel,plot = F, verbose = F, modelNames = "V")
    # clust <- Mclust(y.fit.pred[i,], G = K.sel,plot = F, verbose = F)
    
    # clust <- em(modelName = "V", data = y.fit.pred[i,], parameters = par)
    
    while(is.null(clust)){
      clust <- Mclust(y.fit.pred[i,], G = K.sel,plot = F, verbose = F, modelNames = "E")
      # clust <- em(modelName = "E", data = y.fit.pred[i,], parameters = par)
    }
    
    order = order(clust$parameters$mean)
    mat.order[i,] = order

    mean.temp = clust$parameters$mean[order]
    prob.temp = clust$parameters$pro[order]
    sigma.temp = clust$parameters$variance$sigma[order]

    # mean.temp = clust$parameters$mean
    # prob.temp = clust$parameters$pro
    # sigma.temp = clust$parameters$variance$sigma
    
    if(sum(is.na(sigma.temp)) != 0){
      sigma.temp[order] = sigma.temp[which(!is.na(sigma.temp))]
    }
    
    if(length(sigma.temp) == 1){
      sigma.temp = rep(sigma.temp,K.sel)
    }
    
    for(k in 1:K.sel){
      mat.respons[,k] = prob.temp[k]*dnorm(y.data,mean = mean.temp[k], sd = sqrt(sigma.temp[k]))
    }
    mat.respons = mat.respons/apply(mat.respons, 1, sum)
    
    temp.clust = as.numeric(apply(mat.respons, 1, which.max))
    mat.clust[i,] = temp.clust # important
    # cat(i)
    
    if(km == TRUE){
      kmeans_result <- kmeans(y.fit.pred[i,], centers = K.sel)
      centers = t(kmeans_result$centers)
      order = order(centers[1,])
      centers = centers[,order]
      fixed_means <- data.frame(centers)
      
      mat.clust.k[i,] <- as.numeric(assign_clusters(y.data, fixed_means))
      # cat(i/N)
    }
  }
  
  
  temp = apply(mat.clust,2,function(x) 1 - max(xtabs(~x)/N))
  data.clust = data.frame(dat = y.data, clust.prob = temp)
  
  if(km == TRUE){
    prob.km = apply(mat.clust.k,2,function(x) 1 - max(xtabs(~x)/N))
    data.clust = cbind(data.clust, clust.prob.km = prob.km)
  }
  
  return(data.clust)
  
}

#------------------------------------------------#
# list.dc.possum = temp.DPM.mult
# K.sel = 3
# y.data = y.data.old
# y.data = y.data.2
dc.optim.clust = function(list.dc.possum, y.data, K.sel){

  parm = list.dc.possum[[2]]
  mean = parm$mean[[K.sel]]
  prob = parm$pro[[K.sel]]
  cov.optim = parm$sigmasq[[K.sel]]

  N.samp = nrow(y.data)
  mat.respons = matrix(0,N.samp,K.sel)

  for(k in 1:K.sel){
    mat.respons[,k] = prob[k]*dmvnorm(y.data,mean = mean[,k], sigma = cov.optim[,,k])
  }
  mat.respons = mat.respons/apply(mat.respons, 1, sum)

  optim.clust = apply(mat.respons, 1, which.max)

  data.clust = data.frame(dat = y.data, optim.clust = optim.clust)

  return(data.clust)

}

#-----------------------------------#
#-----------------------------------#

dc.optim.clust.uni = function(list.dc.possum, y.data, K.sel){
  
  parm = list.dc.possum[[2]]
  mean = parm$mean[[K.sel]]
  prob = parm$pro[[K.sel]]
  sigma.optim = parm$sigmasq[[K.sel]]
  
  if(sum(is.na(sigma.optim)) != 0){
    sigma.optim[1:K.sel] = sigma.temp[which(!is.na(sigma.temp))]
  }
  
  N.samp = length(y.data)
  mat.respons = matrix(0,N.samp,K.sel)
  
  for(k in 1:K.sel){
    mat.respons[,k] = prob[k]*dnorm(y.data,mean = mean[k], sd = sigma.optim[k])
  }
  
  mat.respons = mat.respons/apply(mat.respons, 1, sum)
  
  optim.clust = apply(mat.respons, 1, which.max)
  
  data.clust = data.frame(dat = y.data, optim.clust = optim.clust)
  
  return(data.clust)
  
}

# library(mclust)
# library(dplyr)
# y.pred = DPM.acidity[[6]]
# y.data = acidity
# possum_clust = possum_clust.DPM
process_clustering_uni <- function(y.pred, possum_clust, K_star, y.data) {
  # Perform Mclust clustering
  # model <- Mclust(y.pred, G = K_star, modelNames = "E")
  # 
  model <- Mclust(y.pred, G = K_star,plot = F, verbose = F, modelNames = "V")
  # clust <- Mclust(y.fit.pred[i,], G = K.sel,plot = F, verbose = F)
  
  while(is.null(model)){
    model <- Mclust(y.pred, G = K_star,plot = F, verbose = F, modelNames = "E")
  }
  
  pred <- predict(model, newdata = y.data)
  
  # Perform k-means clustering
  kmeans_result <- kmeans(y.pred, centers = K_star)
  
  # Extract and order k-means cluster centers
  centers <- t(kmeans_result$centers)
  center_order <- order(centers[1, ])
  centers <- centers[, center_order]
  
  # Assign clusters based on distances to centers
  distances <- abs(outer(y.data, centers, "-"))
  cluster_assignments <- apply(distances, 1, which.min)
  
  # Combine results into a single data frame
  poss_clust <- cbind(
    possum_clust,
    optim.clust = pred$classification,
    optim.clust.k = cluster_assignments
  )
  
  # Return the processed and arranged data frame
  result <- poss_clust %>% 
    arrange(dat)
  
  return(result)
}


