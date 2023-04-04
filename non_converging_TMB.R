rm(list=ls())
library(ctsmrTMB)
library(ggplot2)
library(ctsmr)
library("plotrix")
setwd("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git")
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
  obj$set_initial_state(init_pars[5], 1*diag(1))
  
  return(obj)}

pars2real = function(x) c(exp(x[1]),x[2],exp(x[3]),exp(x[4])^2)
# testing -----

set.seed(001)
pars = c(theta=10, mu=1.5, sigma_x=1.5, sigma_y=5e-1)

N.sim <- 1000
dt.sim <-1e-3
dt.obs <- 1e-2
x0 <- 3

l <- sim_OU_EM(N.sim, dt.sim, dt.obs, pars, x0)
.data <- l$.data
x <- l$x

plot(.data$y,type='l')

#to easy change the start-parameters c(theta, mu, sigma_x, sigma_y,x0)
init_pars <- c(3, 1, 1e-1, 1e-1, 2)

#to easy change lower bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_lb <- c(1e-5, 3, 1e-10, 1e-10, 1)

#to easy change upper bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_ub <- c(50, 5, 10, 10, 100)

obj <- OU_ctsmrTMB(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub)

nll <- obj$construct_nll(.data, method="ekf")

fitTMB <- obj$estimate(.data, method="ekf", use.hessian=T)

summary(fitTMB)#,extended = TRUE)
pars2real(fitTMB$par.fixed)
fitTMB$sd.fixed



# Check parameter estimates against truth

sum_phillip <- cbind( pars2real(fitTMB$par.fixed), pars )
sum_phillip
# plot one-step predictions, simulated states and observations
t.est = fitTMB$states$mean$prior$t
x.mean = fitTMB$states$mean$prior$x
x.sd = sqrt(fitTMB$states$sd$prior$x)
y = .data$y
ggplot() +
  geom_ribbon(aes(x=t.est, ymin=x.mean-2*x.sd, ymax=x.mean+2*x.sd),fill="grey", alpha=0.9) +
  geom_line(aes(x=t.sim,y=x)) + 
  geom_line(aes(x=t.est, x.mean),col="blue") +
  geom_point(aes(x=t.obs,y=y),col="red",size=2) +
  theme_minimal()


# Check one-step-ahead residuals
plot(fitTMB, use.ggplot=F)

# Use prediction function to get 10-step-ahead state predictions
pred = predict(fitTMB, n.step.ahead=10)

