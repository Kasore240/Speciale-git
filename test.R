#test
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

N.sim <- 1
dt.sim <-1e-3
dt.obs <- 1e-2
x0 <- 3
sx0 <- 2


obj <- OU_ctsmrTMB(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub,sx0)
model <- OU_ctsmr(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub)

  
l <- sim_OU_EM(1000, 1e-3, 1e-2, pars,x0)
.data <- l$.data
  
fitTMB <- obj$estimate(high_sd , method="ekf", use.hessian=T)
fit <- model$estimate(high_sd )

summary(fit)
summary(fitTMB)

plot(.data$t,.data$y,type="l",col="red")
lines(high_sd$t,high_sd$y,col="blue")

high_sd <- .data
