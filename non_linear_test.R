#load librarties
rm(list=ls())
library(ctsmrTMB)
library(ggplot2)
library(ctsmr)

setwd("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git")
####---------Make function to Simulate data using Euler Maruyama sim_OU_EM()-----

sim_OU_sqrt_EM <- function(N.sim, dt.sim, dt.obs, pars,x0){
  
  t.sim = seq(0,1,by=dt.sim)
  dw = rnorm(length(t.sim)-1,sd=sqrt(dt.sim))
  x = x0
  for(i in 1:N.sim) {
    x[i+1] = x[i] + pars[1]*(pars[2]-x[i])*dt.sim *x[i] + pars[3]*dw[i]
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
OU_sqrt_ctsmr <- function(init_pars,init_lb,init_ub){
  
  model <- ctsm$new()
  
  # Add system equations
  model$addSystem( dx ~ theta * (mu-x) * dt *x + exp(log_sigma_x)*dw)
  
  # Add observation equations
  model$addObs(y ~ x)
  
  # Set observation equation variances
  model$setVariance(y ~ (exp(logsigma2_y))^2)
  
  model$setParameter(x0 = c(init=init_pars[5],lb=init_lb[5],ub=init_ub[5]))
  
  # Specify parameter initial values and lower/upper bounds in estimation
  model$setParameter(
    logsigma2_y = log(c(init=init_pars[4], lb=init_lb[4], ub=init_ub[4])),
    log_sigma_x = log(c(init= init_pars[3], lb=init_lb[3], ub=init_ub[3])),
    mu = c(init=init_pars[2], lb=init_lb[2], ub=init_ub[2]),
    theta = c(init = init_pars[1], lb=init_lb[1], ub=init_ub[1])
  )
  
  # set initial value
  
  return(model)
}
#### ---------------- Function to Make new object in ctsmrTMB OU_ctsmrTMB()----
OU_sqrt_ctsmrTMB <- function(init_pars,init_lb,init_ub){
  
  obj = ctsmrTMB$new()
  
 # obj$set_modelname("ornstein_uhlenbeck_sqrt")
  
  # Add system equations
  obj$add_systems(
    dx ~ theta * (mu-x) * dt*x + sigma_x*dw
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
    logsigma_x ~ log(c(init= init_pars[3], lower=init_lb[3], upper=init_ub[3])),
    logsigma_y ~ log(c(init=init_pars[4], lower=init_lb[4], upper=init_ub[4]))
  )
  # Set initial state mean and covariance
  obj$set_initial_state(10, 2*diag(1))
  
  return(obj)}


#### -------------- Change parameters, start guesses and upper/lower bounds ---------
#true parameters 
pars = c(theta=10, mu=1, sigma_x=1.5, sigma_y=1e-2)

#to easy change the start-parameters c(theta, mu, sigma_x, sigma_y,x0)
init_pars <- c(1, 1.5, 1.5,1e-10, 3)

#to easy change lower bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_lb <- c(1e-5, 0, 1e-10, 1e-11, 1)

#to easy change upper bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_ub <- c(50, 5, 20, 10, 100)

N.sim <- 1000
dt.sim <-1e-3
dt.obs <- 1e-2
x0 <- 6

#### -- call functions ---- 

l <- sim_OU_sqrt_EM(N.sim, dt.sim, dt.obs, pars, x0)
.data <- l$.data
.data
x <- l$x
plot(.data$t,.data$y)
model <- OU_sqrt_ctsmr(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub)
obj <- OU_sqrt_ctsmrTMB(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub)


####------------------- Predict in org ctsmr -------------------
# Run the parameter estimation
fit <- model$estimate(.data)
summary(fit,extended = TRUE)
# See the summary of the estimation
summ_org <- summary(fit)
summ_org_exp <- c(fit$xm[5], fit$xm[4], exp(fit$xm[2]),(exp(fit$xm[3])))
summ_org_exp
fit$sd[2:5]



####------------------- Predict in ctsmrTMB -------------------

# Carry out estimation using extended kalman filter method
fitTMB <- obj$estimate(.data, method="ekf", use.hessian=T)

# Check parameter estimates against truth
pars2real = function(x) c(exp(x[1]),x[2],exp(x[3]),exp(x[4]))
summ_TMB <- pars2real(fitTMB$par.fixed)
summary(fitTMB)
### ------ bind results -----------------------

summ_tot <- cbind(summ_TMB,summ_org_exp,pars)
summ_tot


#### ---------------- Plotting -----------------------


t.sim = seq(0,1,by=dt.sim)
t.obs = seq(0,1,by=dt.obs)
### plots org
t.est = seq(0,1,by=dt.obs)
x.mean = tmp$output$pred$y
x.sd = tmp$output$sd$y
par(mfrow=c(2,1))
ggplot() +
  geom_ribbon(aes(x=t.est, ymin=x.mean-2*x.sd, ymax=x.mean+2*x.sd),fill="grey", alpha=0.9) +
  geom_line(aes(x=t.sim,y=x)) + 
  geom_line(aes(x=t.obs, x.mean),col="blue") +
  geom_point(aes(x=t.obs,y=l$.data$y),col="red",size=2) +
  theme_minimal()

### plots TMB
t.est = fitTMB$states$mean$prior$t
x.mean = fitTMB$states$mean$prior$x
x.sd = sqrt(fitTMB$states$sd$prior$x)
ggplot() +
  geom_ribbon(aes(x=t.est, ymin=x.mean-2*x.sd, ymax=x.mean+2*x.sd),fill="grey", alpha=0.9) +
  geom_line(aes(x=t.sim,y=x)) + 
  geom_line(aes(x=t.est, x.mean),col="blue") +
  geom_point(aes(x=t.obs,y=l$.data$y),col="red",size=2) +
  theme_minimal()
plot(fitTMB, use.ggplot=F)

tmp <- predict(fit)[[1]]
str(tmp)

# Calculate the residuals and put them with the data in a data.frame X
.data$residuals <- .data$y - tmp$output$pred$y

# Plot the auto-correlation function and cumulated periodogram in a new window
par(mfrow=c(1,3))
# The blue lines indicates the 95% confidence interval, meaning that if it is
#  white noise, then approximately 19 out of 20 lag correlations will be inside.
acf(.data$residuals, lag.max=8*24)
# The periodogram is the estimated energy spectrum in the signal
spec.pgram(.data$residuals)
# The cumulated periodogram 
cpgram(.data$residuals)
#----------------------------------------------------------------
plot(.data$residuals,type='l')




#### --------------- simulation 200 fits ------------------ 
set.seed(0807)

N.sim<-200
parms_org <- matrix(data=NA,nrow=200,ncol=4)
sd_org <- matrix(data = NA, nrow = 200, ncol = 1)

parms_TMB <- matrix(data=NA,nrow=200,ncol=4)
sd_TMB <- matrix(data = NA, nrow = 200, ncol = 1)

for (i in 1:N.sim){
  l <- sim_OU_EM(N.sim, dt.sim, dt.obs, pars, x0)
  .data <- l$.data
  .data
  x <- l$x
  
  model <- OU_ctsmr(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub)
  obj <- OU_ctsmrTMB(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub)
  
  fit <- model$estimate(.data)
  
  summ_org <- summary(fit)
  summ_org_exp <- c(fit$xm[5], fit$xm[4], exp(fit$xm[2]),(exp(fit$xm[3])))
  parms_org[i,] <- summ_org_exp
  tmp <- predict(fit)[[1]]
  sd_org[i] <- tmp$output$sd$y[101]
  
  
  # Carry out estimation using extended kalman filter method
  fitTMB <- obj$estimate(.data, method="ekf", use.hessian=T)
  
  # Check parameter estimates against truth
  pars2real = function(x) c(exp(x[1]),x[2],exp(x[3]),exp(x[4])^2)
  summ_TMB <- pars2real(fitTMB$par.fixed) 
  parms_TMB[i,] <- summ_TMB
  sd_TMB[i] <- sqrt(fitTMB$states$sd$prior$x)[101]
}
#### --------- dist plots
par(mfrow=c(2,2))
plot(density(parms_TMB[,1]),col="blue",ylim=c(0,0.1), main="theta")
lines(density(parms_org[,1]))

plot(density(parms_TMB[,2]),col="blue",main="mu")
lines(density(parms_org[,2]))

plot(density(parms_TMB[,3]),col="blue", main = "sigma_x")
lines(density(parms_org[,3]))

plot(density(parms_TMB[,4]),col="blue",ylim=c(0,60), main = "sigma_y")
lines(density(parms_org[,4]))

#### -----

summ_org
fitTMB$par.fixed
summ_TMB
