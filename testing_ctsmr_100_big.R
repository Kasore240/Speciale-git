
setwd("/zhome/1e/0/121927/Speciale-git/")
rm(list=ls())
library(ctsmrTMB)
library(ctsmr)
setwd("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git")

## Source functions
sapply(dir("Funktioner",full.names=TRUE), source)
# Parameters -----

N.sim <- 1
dt.sim <-1e-5
dt.obs <- 1e-3
sx0 <- 2

# ---------------------test influence of sigma_x  ------
noise <- seq(1e-1,5,0.49)


parms_org <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))
sd_org <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))

parms_tmb <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))
sd_tmb <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))

set.seed(11)


#to easy change the start-parameters c(theta, mu, sigma_x, sigma_y,x0)
init_pars <- c(1, 1.5, 1e-1, 1e-1, 10)

#to easy change lower bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_lb <- c(1e-5, 0, 1e-10, 1e-10, 1)

#to easy change upper bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_ub <- c(50, 5, 10, 10, 100)

for(j in 1:1){
  pars = c(theta=10, mu=1, sigma_x=noise[j], sigma_y=1e-2)
  
  
  for (i in 1:1){
    
    l <- sim_OU_EM(1000, dt.sim, dt.obs, pars,3)
    .data <- l$.data
    x <- l$x
    
    model <- OU_ctsmr(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub)
    
    fit <- model$estimate(.data)
    
    summ_org <- summary(fit)
    parms_org[i,,j] <- fit$xm[2:5]
    sd_org[i,,j] <- fit$sd[2:5]
    s <-utils::tail(getLoadedDLLs(), 1)
    dyn.unload(s[[1]][["path"]])
    
    
    obj <- OU_ctsmrTMB(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub,sx0)
    
    fitTMB <- obj$estimate(.data, method="ekf", use.hessian=T)
    
    # Check parameter estimates against truth
    summ_TMB <- (fitTMB$par.fixed) 
    parms_tmb[i,,j]  <- summ_TMB
    sd_tmb[i,,j] <- fitTMB$sd.fixed
    
    
  }

  
}

#save(parms_org,sd_org,sd_tmb,parms_tmb,file = "data100big/sigmax_varieret_100.RData")





# --------------------test influence of mu start guess  ------

int_mu <-  seq(0,4.5,0.5)

parms_org <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))
sd_org <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))


parms_tmb <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))
sd_tmb <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))

set.seed(11)

pars = c(theta=10, mu=1, sigma_x=1, sigma_y=1e-2)
#to easy change lower bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_lb <- c(1e-5, -1e-5, 1e-10, 1e-10, 1)

#to easy change upper bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_ub <- c(50, 5, 10, 10, 100)


for(j in 1:10){
  
  #to easy change the start-parameters c(theta, mu, sigma_x, sigma_y,x0)
  init_pars <- c(1, int_mu[j], 1e-1, 1e-1, 3)
  
  for (i in 1:N.sim){
    
    l <- sim_OU_EM(1000, dt.sim, dt.obs, pars, 3)
    .data <- l$.data
    x <- l$x
    
    model <- OU_ctsmr(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub)
    
    fit <- model$estimate(.data)
    
    summ_org <- summary(fit)
    parms_org[i,,j] <- fit$xm[2:5]
    sd_org[i,,j] <- fit$sd[2:5]
    s <-utils::tail(getLoadedDLLs(), 1)
    dyn.unload(s[[1]][["path"]])
    
    obj <- OU_ctsmrTMB(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub,sx0)
    
    fitTMB <- obj$estimate(.data, method="ekf", use.hessian=T)
    
    summ_TMB <- (fitTMB$par.fixed) 
    parms_tmb[i,,j]  <- summ_TMB
    sd_tmb[i,,j] <- fitTMB$sd.fixed
    
  }
  if (j %% 5==0){
    plot(i)
    title('mu start varieret')}
}
save(parms_org,sd_org,sd_tmb,parms_tmb,file = "data100big/Init_mu_vari_100.RData")


# --------------------test influence of theta start guess  ------
int_theta <-  seq(0,18,2)

parms_org <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))
sd_org <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))

parms_tmb <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))
sd_tmb <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))

set.seed(11)

pars = c(theta=10, mu=1, sigma_x=1, sigma_y=1e-2)

#to easy change lower bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_lb <- c(-1e-5, 0, 1e-10, 1e-10, -1)

#to easy change upper bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_ub <- c(50, 5, 10, 10, 100)

for(j in 1:1){
  
  #to easy change the start-parameters c(theta, mu, sigma_x, sigma_y,x0)
  init_pars <- c(1e-10, 1.5, 1e-1, 1e-1, 10)
  
  for (i in 1:N.sim){
    
    l <- sim_OU_EM(1000, dt.sim, dt.obs, pars, 3)
    .data <- l$.data
    x <- l$x
    
    model <- OU_ctsmr(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub)
    
    fit <- model$estimate(.data)
    
    summ_org <- summary(fit)
    parms_org[i,,j] <- fit$xm[2:5]
    
    sd_org[i,,j] <- fit$sd[2:5]
    s <-utils::tail(getLoadedDLLs(), 1)
    dyn.unload(s[[1]][["path"]])
    
    obj <- OU_ctsmrTMB(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub,sx0)
    
    fitTMB <- obj$estimate(.data, method="ekf", use.hessian=T)
    summ_TMB <- (fitTMB$par.fixed) 
    parms_tmb[i,,j]  <- summ_TMB
    sd_tmb[i,,j] <- fitTMB$sd.fixed
    
  }
}
save(parms_org,sd_org,sd_tmb,parms_tmb,file = "data100big/Init_theta_vari_100.RData")


#--------------------test influence of x0 start guess  ------

int_x0 <-  seq(-1,12.5,1.5)

parms_org <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))
sd_org <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))

parms_tmb <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))
sd_tmb <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))

set.seed(11)

pars = c(theta=10, mu=1, sigma_x=1, sigma_y=1e-2)

#to easy change lower bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_lb <- c(1e-5, 0, 1e-10, 1e-10, -1.5)

#to easy change upper bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_ub <- c(50, 5, 10, 10, 100)

for(j in 1:10){
  pars = c(theta=10, mu=1, sigma_x=1, sigma_y=1e-2)
  
  #to easy change the start-parameters c(theta, mu, sigma_x, sigma_y,x0)
  init_pars <- c(1, 1.5, 1e-1, 1e-1, int_x0[j])
  
  for (i in 1:N.sim){
    
    l <- sim_OU_EM(1000, dt.sim, dt.obs, pars, 3)
    .data <- l$.data
    x <- l$x
    
    model <- OU_ctsmr(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub)
    
    fit <- model$estimate(.data)
    
    summ_org <- summary(fit)
    parms_org[i,,j] <- fit$xm[2:5]
    
    sd_org[i,,j] <- fit$sd[2:5]
    s <-utils::tail(getLoadedDLLs(), 1)
    dyn.unload(s[[1]][["path"]])
    
    obj <- OU_ctsmrTMB(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub,sx0)
    
    fitTMB <- obj$estimate(.data, method="ekf", use.hessian=T)
    
    summ_TMB <- (fitTMB$par.fixed) 
    parms_tmb[i,,j]  <- summ_TMB
    sd_tmb[i,,j] <- fitTMB$sd.fixed
    
  }
  if (j %% 5==0){
    plot(i)
    title('x0 varieret')}
}
save(parms_org,sd_org,sd_tmb,parms_tmb,file = "data100big/init_x0_varieret_100.RData")


# -------------------test influence of sig_x start guess ------

int_sig_x <-  seq(1e-10,3,0.3)


parms_org <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))
sd_org <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))
parms_tmb <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))
sd_tmb <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))

set.seed(11)

pars = c(theta=10, mu=1, sigma_x=1, sigma_y=1e-2)

#to easy change lower bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_lb <- c(1e-5, 0, 1e-11, 1e-10, -1.5)

#to easy change upper bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_ub <- c(50, 5, 10, 10, 100)


for(j in 1:10){
  #to easy change the start-parameters c(theta, mu, sigma_x, sigma_y,x0)
  init_pars <- c(1, 1.5, int_sig_x[j], 1e-1, 10)
  
  for (i in 1:N.sim){
    
    l <- sim_OU_EM(1000, dt.sim, dt.obs, pars, 3)
    .data <- l$.data
    x <- l$x
    
    model <- OU_ctsmr(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub)
    
    fit <- model$estimate(.data)
    
    summ_org <- summary(fit)
    parms_org[i,,j] <- fit$xm[2:5]
    
    sd_org[i,,j] <- fit$sd[2:5]
    s <-utils::tail(getLoadedDLLs(), 1)
    dyn.unload(s[[1]][["path"]])
    
    obj <- OU_ctsmrTMB(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub,sx0)
    
    fitTMB <- obj$estimate(.data, method="ekf", use.hessian=T)
    
    summ_TMB <- (fitTMB$par.fixed) 
    parms_tmb[i,,j]  <- summ_TMB
    sd_tmb[i,,j] <- fitTMB$sd.fixed
    
    
  }
  if (j %% 5==0){
    plot(i)
    title('sigma x start varieret')}
}

save(parms_org,sd_org,sd_tmb,parms_tmb,file = "data100big/init_sigmax_varieret_100.RData")

# -------------------test influence of sig_y start guess ------
int_sig_y <-  seq(1e-5,1,1e-1)
parms_org <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))
sd_org <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))

parms_tmb <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))
sd_tmb <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))

set.seed(11)

pars = c(theta=10, mu=1, sigma_x=1, sigma_y=1e-2)

#to easy change lower bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_lb <- c(1e-5, 0, 1e-10, 1e-12, -1.5)

#to easy change upper bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_ub <- c(50, 5, 10, 10, 100)


for(j in 1:10){
  
  #to easy change the start-parameters c(theta, mu, sigma_x, sigma_y,x0)
  init_pars <- c(1, 1.5, 1e-1, int_sig_y[j],10)
  
  for (i in 1:N.sim){
    
    l <- sim_OU_EM(1000, dt.sim, dt.obs, pars, 3)
    .data <- l$.data
    x <- l$x
    
    model <- OU_ctsmr(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub)
    
    fit <- model$estimate(.data)
    
    summ_org <- summary(fit)
    parms_org[i,,j] <- fit$xm[2:5]
    sd_org[i,,j] <- fit$sd[2:5]
    s <-utils::tail(getLoadedDLLs(), 1)
    dyn.unload(s[[1]][["path"]])
    
    obj <- OU_ctsmrTMB(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub,sx0)
    
    fitTMB <- obj$estimate(.data, method="ekf", use.hessian=T)
    
    summ_TMB <- (fitTMB$par.fixed) 
    parms_tmb[i,,j]  <- summ_TMB
    sd_tmb[i,,j] <- fitTMB$sd.fixed
    
  }
  if (j %% 5==0){
    plot(i)
    title('sigma y start varieret')}
}
save(parms_org,sd_org,sd_tmb,parms_tmb,file = "data100big/init_sigmay_varieret_100.RData")