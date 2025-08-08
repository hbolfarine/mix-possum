library(coda)
library(MASS)
library(DPpackage)
library(pracma)

data.bp = read.csv("data_bp.csv")
hist(data.bp$x)
length(data.bp$x)
# data.bp<-data.frame(speeds = galaxies/1000) 

# Initial state
state <- NULL

# MCMC parameters

nburn<-1000
nsave<-1000
nskip<-10
ndisplay<-100
mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)

# aa0=2.01
# Prior
# prior<-list(aa0=5,
#             ab0 = 0.5,
#             kmax = 500,
#             a0=1,
#             b0=1)

#--------------------------------##--------------------------------#
#--------------------------------##--------------------------------#
# Laplace

data.bp = read.csv("examp_laplace.csv")
hist(data.bp$x)
length(data.bp$x)

hist(rgamma(5000,shape = 2, scale = 0.1))
prior<-list(aa0 = 2,
            ab0 = 0.1,
            kmax = 300,
            a0=1,
            b0=1)

set.seed(1820)

# check if the number of components can be changed
fit <- BDPdensity(y=data.bp$x,prior=prior,mcmc=mcmc, ngrid = 500,
                  state=state,status=TRUE)

post.samp.par = fit$save.state$thetasave
write.csv(post.samp.par,row.names = FALSE,"bp_par_lapalce.csv")
n.obs = dim(data.bp)[1]
post.samp.Y = as.matrix(fit$save.state$randsave[,1:n.obs]) 
write.csv(post.samp.Y,row.names = FALSE,"post_samp_Y_laplace.csv")
post.pred = fit$save.state$randsave[,(n.obs+1)]
write.csv(post.pred,row.names = FALSE,"post_pred_laplace.csv")
hist(post.pred,breaks = 100)

table(post.samp.par[,2])
plot(fit$grid,fit$fun)

write.csv(fit$grid,row.names = FALSE,"grid_exemp_laplace.csv")
write.csv(fit$fun,row.names = FALSE,"func_pred_exemp_laplace.csv")

#--------------------------------##--------------------------------#
#--------------------------------##--------------------------------#
# Example 02

# data.bp = read.csv("data_bp_examp_02.csv")
data.bp = read.csv("example_02/data_bp_examp_02_600.csv")

hist(data.bp$x)
length(data.bp$x)
# data.bp<-data.frame(speeds = galaxies/1000) 

# Initial state
state <- NULL

# MCMC parameters

nburn<-1000
nsave<-5000
nskip<-10
ndisplay<-100
mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)

prior<-list(aa0 = 2.01,
            ab0 = 0.1,
            kmax = 1000,
            a0=1,
            b0=1)

set.seed(1820)

# min.y = min(y.data) - 0.1*min(y.data)
# max.y = max(y.data) + 0.1*max(y.data)

# min.y = 10
# max.y = 40

# grid.seq = seq(min.y, max.y, length.out = 500)

# check if the number of components can be changed
fit <- BDPdensity(y=data.bp$x,prior=prior,mcmc=mcmc, ngrid = 500,
                  state=state,status=TRUE)

post.samp.par = fit$save.state$thetasave
write.csv(post.samp.par,row.names = FALSE,"example_02/bp_par_example_02.csv")
n.obs = dim(data.bp)[1]
post.samp.Y = as.matrix(fit$save.state$randsave[,1:n.obs]) 
write.csv(post.samp.Y,row.names = FALSE,"example_02/post_samp_Y_example_02.csv")
post.pred = fit$save.state$randsave[,(n.obs+1)]
write.csv(post.pred,row.names = FALSE,"example_02/post_pred_example_02.csv")
hist(post.pred,breaks = 100)

table(post.samp.par[,2])
plot(fit$grid,fit$fun)
write.csv(fit$fun,row.names = FALSE,"example_02/fun_example_02.csv")

#--------------------------------##--------------------------------#
#--------------------------------##--------------------------------#
# Galaxy data

data.bp = read.csv("bp_galaxy.csv")

hist(data.bp$x)
length(data.bp$x)
# data.bp<-data.frame(speeds = galaxies/1000) 

# Initial state
state <- NULL

# MCMC parameters

nburn<-1000
nsave<-1000
nskip<-10
ndisplay<-100
mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)

prior<-list(aa0 = 2,
            ab0 = 0.01,
            kmax = 500,
            a0=1,
            b0=1)

set.seed(1820)

# check if the number of components can be changed
fit <- BDPdensity(y=data.bp$x,prior=prior,mcmc=mcmc, ngrid = 500,
                  state=state,status=TRUE)

post.samp.par = fit$save.state$thetasave
write.csv(post.samp.par,row.names = FALSE,"galaxy/bp_par_galaxy.csv")
n.obs = dim(data.bp)[1]
post.samp.Y = as.matrix(fit$save.state$randsave[,1:n.obs]) 
write.csv(post.samp.Y,row.names = FALSE,"galaxy/post_samp_Y_galaxy.csv")
post.pred = fit$save.state$randsave[,(n.obs+1)]
write.csv(post.pred,row.names = FALSE,"galaxy/post_pred_galaxy.csv")
hist(post.pred,breaks = 100)

table(post.samp.par[,2])
plot(fit$grid,fit$fun)

write.csv(fit$grid,row.names = FALSE,"galaxy/grid_galaxy.csv")
write.csv(fit$fun,row.names = FALSE,"galaxy/func_galaxy.csv")

#--------------------------------##--------------------------------#
#--------------------------------##--------------------------------#

library(mclust)
data("acidity")

# Initial state
state <- NULL

# MCMC parameters

nburn<-1000
nsave<-5000
nskip<-10
ndisplay<-100
mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)

prior<-list(aa0=2.01,
            ab0=0.01,
            kmax=100,
            a0=1,
            b0=1)

# prior<-list(aa0 = 2,
#             ab0 = 1,
#             kmax = 200,
#             a0=1,
#             b0=1)

set.seed(1820)
n.obs = length(acidity)

# check if the number of components can be changed
fit <- BDPdensity(acidity,prior=prior,mcmc=mcmc, ngrid = 500,
                  state=state,status=TRUE)

post.samp.par = fit$save.state$thetasave
write.csv(post.samp.par,row.names = FALSE,"acidity/bp_par_acidity.csv")
# n.obs = dim(data.bp)[1]
post.samp.Y = as.matrix(fit$save.state$randsave[,1:n.obs]) 
write.csv(post.samp.Y,row.names = FALSE,"acidity/post_samp_Y_acidity.csv")
post.pred = fit$save.state$randsave[,(n.obs+1)]
write.csv(post.pred,row.names = FALSE,"acidity/post_acidity.csv")
hist(post.pred,breaks = 100)

table(post.samp.par[,2])
plot(fit$grid,fit$fun)

write.csv(fit$grid,row.names = FALSE,"acidity/grid_acidity.csv")
write.csv(fit$fun,row.names = FALSE,"acidity/func_acidity.csv")


#--------------------------------##--------------------------------#
#--------------------------------##--------------------------------#


set.seed(1890)
load(file = "enzyme.rda")
n.obs = length(enzyme)

# Initial state
state <- NULL

# MCMC parameters

nburn<-1000
nsave<-5000
nskip<-10
ndisplay<-100
mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)

prior<-list(aa0=1.5,
            ab0=0.01,
            kmax=200,
            a0=1,
            b0=1)

# prior<-list(aa0 = 2,
#             ab0 = 1,
#             kmax = 200,
#             a0=1,
#             b0=1)

# check if the number of components can be changed
fit <- BDPdensity(enzyme,prior=prior,mcmc=mcmc, ngrid = 500,
                  state=state,status=TRUE)

post.samp.par = fit$save.state$thetasave
write.csv(post.samp.par,row.names = FALSE,"enzyme/bp_par_enzyme.csv")
# n.obs = dim(data.bp)[1]
post.samp.Y = as.matrix(fit$save.state$randsave[,1:n.obs]) 
write.csv(post.samp.Y,row.names = FALSE,"enzyme/post_samp_Y_enzyme.csv")
post.pred = fit$save.state$randsave[,(n.obs+1)]
write.csv(post.pred,row.names = FALSE,"enzyme/post_pred_enzyme.csv")
hist(post.pred,breaks = 100)

table(post.samp.par[,2])
plot(fit$grid,fit$fun)

write.csv(fit$grid,row.names = FALSE,"enzyme/grid_enzyme.csv")
write.csv(fit$fun,row.names = FALSE,"enzyme/func_enzyme.csv")

#--------------------------------##--------------------------------#
#--------------------------------##--------------------------------#

# simulation - n = 250

data.bp = read.csv("sim_hellinger/data_bp_01_01.csv")
n.obs = dim(data.bp)[1]

hist(data.bp$x)
length(data.bp$x)

# Initial state
state <- NULL

# MCMC parameters

nburn<-1000
nsave<-1000
nskip<-10
ndisplay<-100
mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)

prior<-list(aa0 = 2,
            ab0 = 1,
            kmax = 500,
            a0=1,
            b0=1)

set.seed(1820)

# check if the number of components can be changed
fit <- BDPdensity(y=data.bp$x,prior=prior,mcmc=mcmc,
                  state=state,status=TRUE)

post.samp.par.250 = fit$save.state$thetasave
write.csv(post.samp.par.250,row.names = FALSE,"sim_hellinger/output_250/bp_par_example_250.csv")
post.samp.Y.250 = as.matrix(fit$save.state$randsave[,1:n.obs]) 
write.csv(post.samp.Y.250,row.names = FALSE,"sim_hellinger/output_250/post_samp_Y_example_250.csv")
post.pred.250 = fit$save.state$randsave[,(n.obs+1)]
write.csv(post.pred.250,row.names = FALSE,"sim_hellinger/output_250/post_pred_example_250.csv")
hist(post.pred,breaks = 100)

table(post.samp.par[,2])
plot(fit$grid,fit$fun)

#--------------------------------##--------------------------------#
#--------------------------------##--------------------------------#
# simulation - n = 500

data.bp = read.csv("sim_hellinger/data_bp_01_02.csv")
n.obs = dim(data.bp)[1]

hist(data.bp$x)
length(data.bp$x)

# Initial state
state <- NULL

# MCMC parameters

nburn<-1000
nsave<-1000
nskip<-10
ndisplay<-100
mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)

prior<-list(aa0 = 2,
            ab0 = 1,
            kmax = 500,
            a0=1,
            b0=1)

set.seed(1820)

# check if the number of components can be changed
fit <- BDPdensity(y=data.bp$x,prior=prior,mcmc=mcmc,
                  state=state,status=TRUE)

post.samp.par.500 = fit$save.state$thetasave
write.csv(post.samp.par.500,row.names = FALSE,"sim_hellinger/output_500/bp_par_example_500.csv")
post.samp.Y.500 = as.matrix(fit$save.state$randsave[,1:n.obs]) 
write.csv(post.samp.Y.500,row.names = FALSE,"sim_hellinger/output_500/post_samp_Y_example_500.csv")
post.pred.500 = fit$save.state$randsave[,(n.obs+1)]
write.csv(post.pred.500,row.names = FALSE,"sim_hellinger/output_500/post_pred_example_500.csv")
hist(post.pred,breaks = 100)

table(post.samp.par[,2])
plot(fit$grid,fit$fun)

#--------------------------------##--------------------------------#
#--------------------------------##--------------------------------#
# simulation - n = 1000

data.bp = read.csv("sim_hellinger/data_bp_01_03.csv")
n.obs = dim(data.bp)[1]

hist(data.bp$x)
length(data.bp$x)

# Initial state
state <- NULL

# MCMC parameters

nburn<-1000
nsave<-1000
nskip<-10
ndisplay<-100
mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)

prior<-list(aa0 = 2,
            ab0 = 3,
            kmax = 500,
            a0=1,
            b0=1)

set.seed(1820)

# check if the number of components can be changed
fit <- BDPdensity(y=data.bp$x,prior=prior,mcmc=mcmc,ngrid= 500,
                  state=state,status=TRUE)

post.samp.par.1000 = fit$save.state$thetasave
write.csv(post.samp.par.1000,row.names = FALSE,"sim_hellinger/output_1000/bp_par_example_1000.csv")
post.samp.Y.1000 = as.matrix(fit$save.state$randsave[,1:n.obs]) 
write.csv(post.samp.Y.1000,row.names = FALSE,"sim_hellinger/output_1000/post_samp_Y_example_1000.csv")
post.pred.1000 = fit$save.state$randsave[,(n.obs+1)]
write.csv(post.pred.1000,row.names = FALSE,"sim_hellinger/output_1000/post_pred_example_1000.csv")
hist(post.pred,breaks = 100)

table(post.samp.par[,2])
plot(fit$grid,fit$fun)

write.csv(fit$grid,row.names = FALSE,"sim_hellinger/output_1000/grid_exemp_02_1000.csv")
write.csv(fit$fun,row.names = FALSE,"sim_hellinger/output_1000/func_exemp_02_1000.csv")


#---------------------#
data.bp

prior2 <- list(a0=2,b0=1,m2=rep(0,1),s2=diag(100000,1),
               psiinv2=solve(diag(0.5,1)),
               nu1=4,nu2=4,tau1=1,tau2=100)

prior1 <- list(alpha=1,m1=rep(0,1),psiinv1=diag(0.5,1),nu1=4,
               tau1=1,tau2=100)

mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)

fit <- DPdensity(y=data.bp$x,prior=prior1,mcmc=mcmc,ngrid= 500,
                 state=state,status=TRUE)

fit$state$alpha

hist(fit$save.state$randsave[,(600+1)])

