### ---- ctsmr ------

#load librarties
rm(list=ls())
library(ctsmrTMB)
library(ggplot2)
library(ctsmr)

setwd("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git")
####---------Make function to Simulate data using Euler Maruyama sim_OU_EM()-----

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
####--------------Make function to make new object in org ctsmr OU_ctsmr() ----
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

#true parameters 
pars = c(theta=10, mu=1, sigma_x=1, sigma_y=1e-2)

#to easy change the start-parameters c(theta, mu, sigma_x, sigma_y,x0)
init_pars <- c(1, 1.5, 1e-1, 1e-1, 10)

#to easy change lower bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_lb <- c(1e-5, 0, 1e-10, 1e-10, 1)

#to easy change upper bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_ub <- c(50, 5, 10, 10, 100)

N.sim <- 1000
dt.sim <-1e-3
dt.obs <- 1e-2
x0 <- 3

#### -- call functions ---- 

l <- sim_OU_EM(N.sim, dt.sim, dt.obs, pars, x0)
.data <- l$.data
.data
x <- l$x
model <- OU_ctsmr(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub)

#### --- Predict ----- 

# Run the parameter estimation
fit <- model$estimate(.data)

# See the summary of the estimation
summ_org <- summary(fit)
summ_org_exp <- c(fit$xm[5], fit$xm[4], exp(fit$xm[2]),(exp(fit$xm[3])))



#### ---- test influence of sigma_x

theta_n <- matrix(data = NA, nrow=10, ncol =2)
mu_n <- matrix(data = NA, nrow=10, ncol =2)
sig_x_n <- matrix(data = NA, nrow=10, ncol =2)
sig_y_n <- matrix(data = NA, nrow=10, ncol =2)
noise <- seq(1e-1,5,0.49)

for(j in 1:10){
pars = c(theta=10, mu=1, sigma_x=noise[j], sigma_y=1e-2)
  
  #to easy change the start-parameters c(theta, mu, sigma_x, sigma_y,x0)
init_pars <- c(1, 1.5, 1e-1, 1e-1, 10)


N.sim<-20
parms_org <- matrix(data=NA,nrow=N.sim,ncol=4)
sd_org <- matrix(data = NA, nrow = N.sim, ncol = 1)


for (i in 1:N.sim){
  l <- sim_OU_EM(1000, dt.sim, dt.obs, pars, x0)
  .data <- l$.data
  x <- l$x
  
  model <- OU_ctsmr(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub)
  
  fit <- model$estimate(.data)
  
  summ_org <- summary(fit)
  summ_org_exp <- c(fit$xm[5], fit$xm[4], exp(fit$xm[2]),(exp(fit$xm[3])))
  parms_org[i,] <- summ_org_exp
  tmp <- predict(fit)[[1]]
  sd_org[i] <- tmp$output$sd$y[101]
}
theta_n[j,1] <- var(parms_org[,1]) 
theta_n[j,2] <- mean(parms_org[,1]) 

mu_n[j,1] <- var(parms_org[,2]) 
mu_n[j,2] <- mean(parms_org[,2]) 

sig_x_n[j,1] <- var(parms_org[,3]) 
sig_x_n[j,2] <- mean(parms_org[,3])

sig_y_n[j,1] <- var(parms_org[,4]) 
sig_y_n[j,2] <- mean(parms_org[,4]) 
}
