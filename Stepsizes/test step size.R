rm(list=ls())
library(ctsmrTMB)
library(ctsmr)
setwd("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git")

## Source functions
sapply(dir("Funktioner",full.names=TRUE), source)

pars = c(theta=10, mu=1, sigma_x=1, sigma_y=1e-2)

#to easy change the start-parameters c(theta, mu, sigma_x, sigma_y,x0)
init_pars <- c(1, 1.5, 1e-1, 1e-1, 10)

#to easy change lower bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_lb <- c(1e-5, 0, 1e-10, 1e-10, 1)

#to easy change upper bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_ub <- c(50, 5, 10, 10, 100)


N.sim <- 1000
dt.sim <-1e-2
dt.obs <- 1e-1
x0 <- 3
sx0 <- 2
yall <- matrix(data=NA,nrow=101,ncol=1000)

parms_tmb <- matrix(data = NA, nrow = 100, ncol = 4)
sd_tmb <- matrix(data = NA, nrow = 100, ncol = 4)
parms_tmb_ts <- matrix(data = NA, nrow = 100, ncol = 4)
sd_tmb_ts <- matrix(data = NA, nrow = 100, ncol = 4)

parms_org <- matrix(data = NA, nrow = 100, ncol = 4)
sd_org <- matrix(data = NA, nrow = 100, ncol = 4)
parms_org_ts <- matrix(data = NA, nrow = 100, ncol = 4)
sd_org_ts <- matrix(data = NA, nrow = 100, ncol = 4)

for(i in 1:10){
  obj <- OU_ctsmrTMB(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub,sx0)
  model <- OU_ctsmr(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub)
  
  l <- sim_OU_EM(1000, 1e-3, 1e-2, pars,x0)
  .data <- l$.data

    fitTMB <- obj$estimate(.data, method="ekf", use.hessian=T)
  summ_TMB <- (fitTMB$par.fixed) 
  parms_tmb[i,]  <- summ_TMB
  sd_tmb[i,] <- fitTMB$sd.fixed
  
  
  
  fit <- model$estimate(.data)
  
  summ_org <- summary(fit)
  parms_org[i,] <- fit$xm[2:5]
  sd_org[i,] <- fit$sd[2:5]
  s <-utils::tail(getLoadedDLLs(), 1)
  dyn.unload(s[[1]][["path"]])
  
  
  l <- sim_OU_EM(1000, dt.sim, dt.obs, pars,x0)
  .data <- l$.data
  
  fitTMB_ts <- obj$estimate(.data, method="ekf", use.hessian=T)
  
  summ_TMB_ts <- (fitTMB_ts$par.fixed) 
  parms_tmb_ts[i,]  <- summ_TMB_ts
  sd_tmb_ts[i,] <- fitTMB_ts$sd.fixed
  
  fit_ts <- model$estimate(.data)
  
  summ_org_ts <- summary(fit_ts)
  parms_org_ts[i,] <- fit_ts$xm[2:5]
  sd_org_ts[i,] <- fit_ts$sd[2:5]
  s <-utils::tail(getLoadedDLLs(), 1)
  dyn.unload(s[[1]][["path"]])
}


par(mfrow=c(2,2))
#theta
hist(exp(parms_tmb[,1]),main ="TMB - big steps", xlab=expression(paste("Estimated ",theta)))
hist(exp(parms_tmb_ts[,1]),main ="TMB - small steps", xlab=expression(paste("Estimated ",theta)))
hist(parms_org[,4],main = "org - big steps", xlab=expression(paste("Estimated ",theta)))
hist(parms_org_ts[,4],main = "org - small steps", xlab=expression(paste("Estimated ",theta)))
mtext(expression(paste('Distribution of ', theta)),side=3,outer=TRUE, line=-1.5,cex.main=1.4)

par(mfrow=c(2,2))
#mu
hist((parms_tmb[,2]),main ="TMB - big steps", xlab=expression(paste("Estimated ",mu)))
hist((parms_tmb_ts[,2]),main ="TMB - small steps", xlab=expression(paste("Estimated ",mu)))
hist(parms_org[,3],main = "org - big steps", xlab=expression(paste("Estimated ",mu)))
hist(parms_org_ts[,3],main = "org - small steps", xlab=expression(paste("Estimated ",mu)))
mtext(expression(paste('Distribution of ', mu)),side=3,outer=TRUE, line=-1.5,cex.main=1.4)

par(mfrow=c(2,2))
#sigma x
hist(exp(parms_tmb[,3]),main ="TMB - big steps", xlab=expression(paste("Estimated ",sigma[x])))
hist(exp(parms_tmb_ts[,3]),main ="TMB - small steps", xlab=expression(paste("Estimated ",sigma[x])))
hist(exp(parms_org[,1]),main = "org - big steps", xlab=expression(paste("Estimated ",sigma[x])))
hist(exp(parms_org_ts[,1]),main = "org - small steps", xlab=expression(paste("Estimated ",sigma[x])))
mtext(expression(paste('Distribution of ', sigma[x])),side=3,outer=TRUE, line=-1.5,cex.main=1.4)

par(mfrow=c(2,2))
#sigma y
hist(exp(parms_tmb[,4]),main ="TMB - big steps", xlab=expression(paste("Estimated ",sigma[y])))
hist(exp(parms_tmb_ts[,4]),main ="TMB - small steps", xlab=expression(paste("Estimated ",sigma[y])))
hist(exp(parms_org[,2]),main = "org - big steps", xlab=expression(paste("Estimated ",sigma[y])))
hist(exp(parms_org_ts[,2]),main = "org - small steps", xlab=expression(paste("Estimated ",sigma[y])))
mtext(expression(paste('Distribution of ', sigma[y])),side=3,outer=TRUE, line=-1.5,cex.main=1.4)

#save(parms_org,parms_org_ts,sd_org,sd_org_ts,parms_tmb,parms_tmb_ts,sd_tmb,sd_tmb_ts, file = 'Stepsizes/stepsize_100sim.Rdata' )

it <-
i <- 2
p <- 1
(c(sd(exp(parms_tmb[,i])),mean(exp(parms_tmb[,i])),sd(exp(parms_tmb_ts[,i])),mean(exp(parms_tmb_ts[,i]))))
#coverage 

    low =(c(parms_tmb_ts[,i]) - (2*c(sd_tmb_ts[,i])))
    up =(c(parms_tmb_ts[,i]) + (2*c(sd_tmb_ts[,i])))
    x <- sum(p< up & p> low,na.rm=T)
    b <- binom.test(x,100,0.95,conf.level = 0.95)


