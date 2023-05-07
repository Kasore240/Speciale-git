rm(list=ls())
library(ctsmrTMB)
library(ctsmr)
setwd("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git")

## Source functions
sapply(dir("Funktioner",full.names=TRUE), source)

#### -------------- Change parameters, start guesses and upper/lower bounds ---------
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
sim <- 20
po <- matrix(data = NA, nrow=sim, ncol = 4)
pt <- matrix(data = NA, nrow=sim, ncol = 4)

for (i in 1:sim){

l <- sim_OU_EM(1000, dt.sim, dt.obs, pars,x0)
.data <- l$.data
x <- l$x
model <- OU_ctsmr(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub)
obj <- OU_ctsmrTMB(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub,2)

fit <- model$estimate(.data)
fitTMB <- obj$estimate(.data, method="ekf", use.hessian=T,ode.timestep = 1e-4)
po[i,] <-fit$xm[2:5]
pt[i,] <- fitTMB$par.fixed
}
save(po,pt,file ='fordeling.Rdata')
