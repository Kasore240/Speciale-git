
#load librarties
rm(list=ls())
library(ctsmrTMB)
library(ggplot2)
library(ctsmr)
library("plotrix")
setwd("/zhome/1e/0/121927/Speciale-git/")
#setwd("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git")
# Make function to Simulate data using Euler Maruyama sim_OU_EM() =========

sim_OU_EM <- function(N.sim, dt.sim, dt.obs, pars,x0){
  
  t.sim = seq(0,1,by=dt.sim)
  dw = rnorm(length(t.sim)-1,sd=sqrt(dt.sim))
  x = x0
  for(i in 1:N.sim) {
    x[i+1] = x[i] + pars[1]*(pars[2]-x[i])*dt.sim + pars[3]*dw[i]
  }
  
  # Extract observations and add noise
  t.obs = seq(0,1,by=dt.obs)
  y = x[t.sim %in% t.obs] + pars[4] * rnorm(length(t.obs))
  
  
  # Create data
  .data = data.frame(
    t = t.obs,
    y = y
  )
  l1 <- list(.data=.data,x=x)
  return(l1)}
# Make function to make new object in org ctsmr OU_ctsmr() ----
OU_ctsmr <- function(init_pars,init_lb,init_ub){
  
  model <- ctsm$new()
  
  # Add system equations
  model$addSystem( dx ~ theta * (mu-x) * dt + exp(log_sigma_x)*dw)
  
  # Add observation equations
  model$addObs(y ~ x)
  
  # Set observation equation variances
  model$setVariance(y ~ (exp(logsigma2_y))^2)
  
  # Specify parameter initial values and lower/upper bounds in estimation
  model$setParameter(
    theta = c(init = init_pars[1], lb=init_lb[1], ub=init_ub[1]),
    mu = c(init=init_pars[2], lb=init_lb[2], ub=init_ub[2]),
    log_sigma_x = log(c(init= init_pars[3], lb=init_lb[3], ub=init_ub[3])),
    logsigma2_y = log(c(init=init_pars[4], lb=init_lb[4], ub=init_ub[4]))
  )
  
  # set initial value
  model$setParameter(x0 = c(init=init_pars[5],lb=init_lb[5],ub=init_ub[5]))
  return(model)
}


# Parameters -----

N.sim <- 100
dt.sim <-1e-3
dt.obs <- 1e-2


# ---------------------test influence of sigma_x  ------
noise <- seq(1e-1,5,0.49)

theta_n <- matrix(data = NA, nrow=10, ncol =2)
mu_n <- matrix(data = NA, nrow=10, ncol =2)
sig_x_n <- matrix(data = NA, nrow=10, ncol =2)
sig_y_n <- matrix(data = NA, nrow=10, ncol =2)

sd_mean <- matrix(data = NA, nrow = 10, ncol = 4)
parms_org <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))
sd_org <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))

set.seed(11)


#to easy change the start-parameters c(theta, mu, sigma_x, sigma_y,x0)
init_pars <- c(1, 1.5, 1e-1, 1e-1, 10)

#to easy change lower bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_lb <- c(1e-5, 0, 1e-10, 1e-10, 1)

#to easy change upper bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_ub <- c(50, 5, 10, 10, 100)

for(j in 1:10){
  pars = c(theta=10, mu=1, sigma_x=noise[j], sigma_y=1e-2)

  
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

  }
  theta_n[j,1] <- var(parms_org[,1,j]) 
  theta_n[j,2] <- mean(parms_org[,1,j]) 
  
  mu_n[j,1] <- var(parms_org[,2,j]) 
  mu_n[j,2] <- mean(parms_org[,2,j]) 
  
  sig_x_n[j,1] <- var(parms_org[,3,j]) 
  sig_x_n[j,2] <- mean(parms_org[,3,j])
  
  sig_y_n[j,1] <- var(parms_org[,4,j]) 
  sig_y_n[j,2] <- mean(parms_org[,4,j]) 
  
  sd_mean[j,1] <- mean(sd_org[,1,j])
  sd_mean[j,2] <- mean(sd_org[,2,j])
  sd_mean[j,3] <- mean(sd_org[,3,j])
  sd_mean[j,4] <- mean(sd_org[,4,j])
}
a_org <- cbind(theta_n,mu_n,sig_x_n,sig_y_n)
save(a_org,parms_org,sd_org,sd_mean,file = "data100/sigmax_varieret_100.RData")


# ---------------------test influence of sigma_y  ------
noise <-  seq(1e-3,4.5,0.49)

theta_n <- matrix(data = NA, nrow=10, ncol =2)
mu_n <- matrix(data = NA, nrow=10, ncol =2)
sig_x_n <- matrix(data = NA, nrow=10, ncol =2)
sig_y_n <- matrix(data = NA, nrow=10, ncol =2)

sd_mean <- matrix(data = NA, nrow = 10, ncol = 4)
parms_org <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))
sd_org <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))

set.seed(11)

#to easy change the start-parameters c(theta, mu, sigma_x, sigma_y,x0)
init_pars <- c(1, 1.5, 1e-1, 1e-1, 10)

#to easy change lower bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_lb <- c(1e-5, 0, 1e-10, 1e-10, 1)

#to easy change upper bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_ub <- c(50, 5, 10, 10, 100)

for(j in 1:10){
  pars = c(theta=10, mu=1, sigma_x=1, sigma_y=noise[j])
  
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
    
  }
  theta_n[j,1] <- var(parms_org[,1,j]) 
  theta_n[j,2] <- mean(parms_org[,1,j]) 
  
  mu_n[j,1] <- var(parms_org[,2,j]) 
  mu_n[j,2] <- mean(parms_org[,2,j]) 
  
  sig_x_n[j,1] <- var(parms_org[,3,j]) 
  sig_x_n[j,2] <- mean(parms_org[,3,j])
  
  sig_y_n[j,1] <- var(parms_org[,4,j]) 
  sig_y_n[j,2] <- mean(parms_org[,4,j]) 
  
  sd_mean[j,1] <- mean(sd_org[,1,j])
  sd_mean[j,2] <- mean(sd_org[,2,j])
  sd_mean[j,3] <- mean(sd_org[,3,j])
  sd_mean[j,4] <- mean(sd_org[,4,j])
}
a_org <- cbind(theta_n,mu_n,sig_x_n,sig_y_n)
save(a_org,parms_org,sd_org,sd_mean,file = "data100/sigmay_varieret_100.RData")


# --------------------test influence of mu start guess  ------

int_mu <-  seq(0,4.5,0.5)

theta_n <- matrix(data = NA, nrow=10, ncol =2)
mu_n <- matrix(data = NA, nrow=10, ncol =2)
sig_x_n <- matrix(data = NA, nrow=10, ncol =2)
sig_y_n <- matrix(data = NA, nrow=10, ncol =2)

sd_mean <- matrix(data = NA, nrow = 10, ncol = 4)
parms_org <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))
sd_org <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))

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
    
  }
  theta_n[j,1] <- var(parms_org[,1,j]) 
  theta_n[j,2] <- mean(parms_org[,1,j]) 
  
  mu_n[j,1] <- var(parms_org[,2,j]) 
  mu_n[j,2] <- mean(parms_org[,2,j]) 
  
  sig_x_n[j,1] <- var(parms_org[,3,j]) 
  sig_x_n[j,2] <- mean(parms_org[,3,j])
  
  sig_y_n[j,1] <- var(parms_org[,4,j]) 
  sig_y_n[j,2] <- mean(parms_org[,4,j]) 
  
  sd_mean[j,1] <- mean(sd_org[,1,j])
  sd_mean[j,2] <- mean(sd_org[,2,j])
  sd_mean[j,3] <- mean(sd_org[,3,j])
  sd_mean[j,4] <- mean(sd_org[,4,j])
}
a_org <- cbind(theta_n,mu_n,sig_x_n,sig_y_n)
save(a_org,parms_org,sd_org,sd_mean,file = "data100/Init_mu_vari_100.RData")


# --------------------test influence of theta start guess  ------
int_theta <-  seq(0,18,2)

theta_n <- matrix(data = NA, nrow=10, ncol =2)
mu_n <- matrix(data = NA, nrow=10, ncol =2)
sig_x_n <- matrix(data = NA, nrow=10, ncol =2)
sig_y_n <- matrix(data = NA, nrow=10, ncol =2)

sd_mean <- matrix(data = NA, nrow = 10, ncol = 4)
parms_org <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))
sd_org <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))

set.seed(11)

pars = c(theta=10, mu=1, sigma_x=1, sigma_y=1e-2)

#to easy change lower bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_lb <- c(-1e-5, 0, 1e-10, 1e-10, -1)

#to easy change upper bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_ub <- c(50, 5, 10, 10, 100)

for(j in 1:10){
  
  #to easy change the start-parameters c(theta, mu, sigma_x, sigma_y,x0)
  init_pars <- c(int_theta[j], 1.5, 1e-1, 1e-1, 10)
  
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
    
  }
  theta_n[j,1] <- var(parms_org[,1,j]) 
  theta_n[j,2] <- mean(parms_org[,1,j]) 
  
  mu_n[j,1] <- var(parms_org[,2,j]) 
  mu_n[j,2] <- mean(parms_org[,2,j]) 
  
  sig_x_n[j,1] <- var(parms_org[,3,j]) 
  sig_x_n[j,2] <- mean(parms_org[,3,j])
  
  sig_y_n[j,1] <- var(parms_org[,4,j]) 
  sig_y_n[j,2] <- mean(parms_org[,4,j]) 
  
  sd_mean[j,1] <- mean(sd_org[,1,j])
  sd_mean[j,2] <- mean(sd_org[,2,j])
  sd_mean[j,3] <- mean(sd_org[,3,j])
  sd_mean[j,4] <- mean(sd_org[,4,j])
}
a_org <- cbind(theta_n,mu_n,sig_x_n,sig_y_n)
save(a_org,parms_org,sd_org,sd_mean,file = "data100/Init_theta_vari_100.RData")


#--------------------test influence of x0 start guess  ------

int_x0 <-  seq(-1,12.5,1.5)

theta_n <- matrix(data = NA, nrow=10, ncol =2)
mu_n <- matrix(data = NA, nrow=10, ncol =2)
sig_x_n <- matrix(data = NA, nrow=10, ncol =2)
sig_y_n <- matrix(data = NA, nrow=10, ncol =2)

sd_mean <- matrix(data = NA, nrow = 10, ncol = 4)
parms_org <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))
sd_org <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))

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
    
  }
  theta_n[j,1] <- var(parms_org[,1,j]) 
  theta_n[j,2] <- mean(parms_org[,1,j]) 
  
  mu_n[j,1] <- var(parms_org[,2,j]) 
  mu_n[j,2] <- mean(parms_org[,2,j]) 
  
  sig_x_n[j,1] <- var(parms_org[,3,j]) 
  sig_x_n[j,2] <- mean(parms_org[,3,j])
  
  sig_y_n[j,1] <- var(parms_org[,4,j]) 
  sig_y_n[j,2] <- mean(parms_org[,4,j]) 
  
  sd_mean[j,1] <- mean(sd_org[,1,j])
  sd_mean[j,2] <- mean(sd_org[,2,j])
  sd_mean[j,3] <- mean(sd_org[,3,j])
  sd_mean[j,4] <- mean(sd_org[,4,j])
}
a_org <- cbind(theta_n,mu_n,sig_x_n,sig_y_n)
save(a_org,parms_org,sd_org,sd_mean,file = "data100/init_x0_varieret_100.RData")


# -------------------test influence of sig_x start guess ------

int_sig_x <-  seq(1e-10,3,0.3)

theta_n <- matrix(data = NA, nrow=10, ncol =2)
mu_n <- matrix(data = NA, nrow=10, ncol =2)
sig_x_n <- matrix(data = NA, nrow=10, ncol =2)
sig_y_n <- matrix(data = NA, nrow=10, ncol =2)

sd_mean <- matrix(data = NA, nrow = 10, ncol = 4)
parms_org <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))
sd_org <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))

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
    
  }
  theta_n[j,1] <- var(parms_org[,1,j]) 
  theta_n[j,2] <- mean(parms_org[,1,j]) 
  
  mu_n[j,1] <- var(parms_org[,2,j]) 
  mu_n[j,2] <- mean(parms_org[,2,j]) 
  
  sig_x_n[j,1] <- var(parms_org[,3,j]) 
  sig_x_n[j,2] <- mean(parms_org[,3,j])
  
  sig_y_n[j,1] <- var(parms_org[,4,j]) 
  sig_y_n[j,2] <- mean(parms_org[,4,j]) 
  
  sd_mean[j,1] <- mean(sd_org[,1,j])
  sd_mean[j,2] <- mean(sd_org[,2,j])
  sd_mean[j,3] <- mean(sd_org[,3,j])
  sd_mean[j,4] <- mean(sd_org[,4,j])
}
a_org <- cbind(theta_n,mu_n,sig_x_n,sig_y_n)

save(a_org,parms_org,sd_org,sd_mean,file = "data100/init_sigmax_varieret_100.RData")

# -------------------test influence of sig_y start guess ------
int_sig_y <-  seq(1e-5,1,1e-1)

theta_n <- matrix(data = NA, nrow=10, ncol =2)
mu_n <- matrix(data = NA, nrow=10, ncol =2)
sig_x_n <- matrix(data = NA, nrow=10, ncol =2)
sig_y_n <- matrix(data = NA, nrow=10, ncol =2)

sd_mean <- matrix(data = NA, nrow = 10, ncol = 4)
parms_org <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))
sd_org <- array(rep(NaN, N.sim*4*10), dim=c(N.sim, 4, 10))

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
    
  }
  theta_n[j,1] <- var(parms_org[,1,j]) 
  theta_n[j,2] <- mean(parms_org[,1,j]) 
  
  mu_n[j,1] <- var(parms_org[,2,j]) 
  mu_n[j,2] <- mean(parms_org[,2,j]) 
  
  sig_x_n[j,1] <- var(parms_org[,3,j]) 
  sig_x_n[j,2] <- mean(parms_org[,3,j])
  
  sig_y_n[j,1] <- var(parms_org[,4,j]) 
  sig_y_n[j,2] <- mean(parms_org[,4,j]) 
  
  sd_mean[j,1] <- mean(sd_org[,1,j])
  sd_mean[j,2] <- mean(sd_org[,2,j])
  sd_mean[j,3] <- mean(sd_org[,3,j])
  sd_mean[j,4] <- mean(sd_org[,4,j])
}
a_org <- cbind(theta_n,mu_n,sig_x_n,sig_y_n)
save(a_org,parms_org,sd_org,sd_mean,file = "data100/init_sigmay_varieret_100.RData")
