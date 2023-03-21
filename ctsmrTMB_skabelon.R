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
#### ---------------- Function to Make new object in ctsmrTMB OU_ctsmrTMB()----
OU_ctsmrTMB <- function(init_pars,init_lb,init_ub){
  
  obj = ctsmrTMB$new()
  
  obj$set_modelname("ornstein_uhlenbeck")
  
  # Add system equations
  obj$add_systems(
    dx ~ theta * (mu-x) * dt + sigma_x*dw
  )
  
  # Add observation equations
  obj$add_observations(
    y ~ x
  )
  
  # Set observation equation variances
  obj$add_observation_variances(y ~ sigma_y^2)
  
  # Specify algebraic relations   -- er de her mon overflødige?
  obj$add_algebraics(
    theta ~ exp(logtheta),
    sigma_x ~ exp(logsigma_x),
    sigma_y ~ exp(logsigma_y)
  )
  
  # Specify parameter initial values and lower/upper bounds in estimation
  obj$add_parameters(
    logtheta ~ log(c(init = init_pars[1], lower=init_lb[1], upper=init_ub[1])),
    mu ~ c(init=init_pars[2], lower=init_lb[1], upper=init_ub[2]),
    logsigma_x ~ log(c(init= init_pars[3], lower=init_lb[1], upper=init_ub[3])),
    logsigma_y ~ log(c(init=init_pars[4], lower=init_lb[1], upper=init_ub[4]))
  )
  # Set initial state mean and covariance
  obj$set_initial_state(3.5, 1e-1*diag(1))
  
  return(obj)}
###### ----- init param ---------------------
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
parms_TMB <- matrix(data=NA,nrow=N.sim,ncol=4)
sd_TMB <- matrix(data = NA, nrow = N.sim, ncol = 1)


for (i in 1:N.sim){
  l <- sim_OU_EM(1000, dt.sim, dt.obs, pars, x0)
  .data <- l$.data
  x <- l$x
  
  obj <- OU_ctsmrTMB(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub)
  fitTMB <- obj$estimate(.data, method="ekf", use.hessian=T)
  
  # Check parameter estimates against truth
  pars2real = function(x) c(exp(x[1]),x[2],exp(x[3]),exp(x[4])^2)
  summ_TMB <- pars2real(fitTMB$par.fixed) 
  parms_TMB[i,] <- summ_TMB
  sd_TMB[i] <- sqrt(fitTMB$states$sd$prior$x)[101]
}
theta_n[j,1] <- var(parms_TMB[,1]) 
theta_n[j,2] <- mean(parms_TMB[,1]) 

mu_n[j,1] <- var(parms_TMB[,2]) 
mu_n[j,2] <- mean(parms_TMB[,2]) 

sig_x_n[j,1] <- var(parms_TMB[,3]) 
sig_x_n[j,2] <- mean(parms_TMB[,3])

sig_y_n[j,1] <- var(parms_TMB[,4]) 
sig_y_n[j,2] <- mean(parms_TMB[,4]) 
}
