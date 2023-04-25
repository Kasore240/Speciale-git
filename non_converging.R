#non converging 
#load librarties
rm(list=ls())
library(ctsmrTMB)
library(ggplot2)
library(ctsmr)
library("plotrix")
setwd("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git")
# Make function to Simulate data using Euler Maruyama sim_OU_EM() =========

sim_OU_EM <- function(N.sim, dt.sim, dt.obs, pars,x0){
  # pars = c(theta=10, mu=1, sigma_x=1, sigma_y=2e-1)
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
set.seed(001)
pars = c(theta=10, mu=1, sigma_x=1, sigma_y=2e-1)

N.sim <- 1000
dt.sim <-1e-3
dt.obs <- 1e-2
x0 <- 3

l <- sim_OU_EM(N.sim, dt.sim, dt.obs, pars, x0)
.data <- l$.data
x <- l$x

plot(.data$y,type='l')

#to easy change the start-parameters c(theta, mu, sigma_x, sigma_y,x0)
init_pars <- c(3, 3, 1, 1, 10)

#to easy change lower bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_lb <- c(1, 2, 1e-10, 1e-10, 1)

#to easy change upper bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_ub <- c(15, 5, 5, 5, 15)


model2 <- OU_ctsmr(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub)

# Run the parameter estimation
fit <- model2$estimate(.data)
fit$message


# See the summary of the estimation
summary(fit, extended = TRUE)

sum_ctsr <- summary(fit, extended = TRUE)
sum_ctsr$coefficients[2:3,'Estimate'] <- exp(sum_ctsr$coefficients[2:3,'Estimate'])
sum_ctsr$coefficients[3,'Estimate'] <- sqrt(sum_ctsr$coefficients[3,'Estimate'])

#----------------------------------------------------------------
sum_ctsr$coefficients[,'Estimate'] 

tmp <- predict(fit)[[1]]
str(tmp)

# Calculate the residuals and put them with the data in a data.frame X
.data$residuals <- .data$y - tmp$output$pred$y

# Plot the auto-correlation function and cumulated periodogram in a new window

# The blue lines indicates the 95% confidence interval, meaning that if it is
#  white noise, then approximately 19 out of 20 lag correlations will be inside.
acf(.data$residuals, lag.max=8*24)
# The periodogram is the estimated energy spectrum in the signal
spec.pgram(.data$residuals)
# The cumulated periodogram 
cpgram(.data$residuals)
#----------------------------------------------------------------
plot(.data$residuals,type='l')

parsfrompr = c(theta=sum_ctsr$coefficients[5,'Estimate'], mu=sum_ctsr$coefficients[4,'Estimate'], sigma_x=sum_ctsr$coefficients[2,'Estimate'], sigma_y=sum_ctsr$coefficients[3,'Estimate'])

lfp <- sim_OU_EM(N.sim, dt.sim, dt.obs, parsfrompr, sum_ctsr$coefficients[1,'Estimate'])
.datafp <- lfp$.data
xfp <- lfp$x

plot(tmp$output$pred$y,type='l',col='red') 
lines(1:101,.data$y,col='blue')
lines(1:101,.datafp$y,col='black')
