#load librarties
rm(list=ls())
library(ctsmrTMB)
library(ggplot2)
library(ctsmr)

setwd("C:/Users/bruger/OneDrive/Skrivebord/Speciale")
####---------Make function to Simulate data using Euler Maruyama sim_OU_EM()-----
set.seed(10)
sim_OU_EM <- function(N.sim, dt.sim,dt.obs, pars,x0){

  t.sim = seq(0,1,by=dt.sim)
  dw = rnorm(length(t.sim)-1,sd=sqrt(dt.sim))
  x = x0
  for(i in 1:N.sim) {
    x[i+1] = x[i] + pars[1]*(pars[2]-x[i])*dt.sim + pars[3]*dw[i]
  }
  
  # Extract observations and add noise
  t.obs = seq(0,1,by=dt.obs)
  y = x[t.sim %in% t.obs] + sqrt(pars[4]) * rnorm(length(t.obs))
  
  
  # Create data
  .data = data.frame(
    t = t.obs,
    y = y
  )
  return(.data)}
####--------------Make function to make new object in org ctsmr OU_ctsmr() ----
OU_ctsmr <- function(init_pars,init_lb,init_ub){
  
  model <- ctsm$new()
  
  # Add system equations
  model$addSystem(dx ~ theta * (mu-x) * dt + exp(log_sigma_x)*dw)
  
  # Add observation equations
  model$addObs(y ~ x)
  
  # Set observation equation variances
  model$setVariance(y ~ exp(logsigma2_y))
  
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
#### ---------------- Function to Make new object in ctsmrTMB OU_ctsmrTMB()----
OU_ctsmrTMB <- function(init_pars,init_lb,init_ub){
  
  obj = ctsmrTMB$new()
  
  obj$set_modelname("ornstein_uhlenbeck")
  
  # Add system equations
  obj$add_systems(dx ~ exp(logtheta) * (mu-x) * dt + exp(logsigma_x)*dw)
  
  # Add observation equations
  obj$add_observations(y ~ x)
  
  # Set observation equation variances
  obj$add_observation_variances(y ~ exp(logsigma_y))
  
  # Specify algebraic relations   -- er de her mon overflødige?
  #obj$add_algebraics(
  # theta ~ exp(logtheta),
  #  sigma_x ~ exp(logsigma_x),
  # sigma_y ~ exp(logsigma_y) )init_lb[1], ub=init_ub[1]
  
  # Specify parameter initial values and lower/upper bounds in estimation
  obj$add_parameters(
    logtheta ~ log(c(init = init_pars[1], lower=init_lb[1], upper=init_ub[1])),
    mu ~ c(init=init_pars[2], lower=init_lb[1], upper=init_ub[2]),
    logsigma_x ~ log(c(init= init_pars[3], lower=init_lb[1], upper=init_ub[3])),
    logsigma_y ~ log(c(init=init_pars[4], lower=init_lb[1], upper=init_ub[4]))
  )
  # Set initial state mean and covariance
  obj$set_initial_state(init_pars[5], 1e-2*diag(1))
  
  return(obj)}


#### -------------- Change parameters, start guesses and upper/lower bounds ---------
#true parameters 
pars = c(theta=10, mu=1, sigma_x=1, sigma_y=1e-2)

#to easy change the start-parameters c(theta, mu, sigma_x, sigma_y,x0)
init_pars <- c(1, 1.5, 1e-1, 1e-1, 10)

#to easy change lower bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_lb <- c(1e-5, 0, 1e-10, 1e-10, 5)

#to easy change upper bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_ub <- c(50, 5, 10, 10, 100)

#### -- call functions ---- 

.data <- sim_OU_EM(1000,1e-3,1e-2,pars,3)
model <- OU_ctsmr(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub)
obj <- OU_ctsmrTMB(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub)


####------------------- Predict in org ctsmr -------------------
# Run the parameter estimation
fit <- model$estimate(.data)

# See the summary of the estimation
summ_org <- summary(fit)
summ_org_exp <- c(fit$xm[5], fit$xm[4], exp(fit$xm[2]),exp(fit$xm[3]))

####------------------- Predict in ctsmrTMB -------------------

# Carry out estimation using extended kalman filter method
fitTMB <- obj$estimate(.data, method="ekf", use.hessian=T)

# Check parameter estimates against truth
pars2real = function(x) c(exp(x[1]),x[2],exp(x[3]),exp(x[4]))
summ_TMB <- pars2real(fitTMB$par.fixed)

### ------ bind results -----------------------

summ_tot <- cbind(summ_TMB,summ_org_exp,pars)
summ_tot

