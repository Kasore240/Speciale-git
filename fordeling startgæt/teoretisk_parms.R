rm(list=ls())
library(ctsmrTMB)
library(ctsmr)
setwd("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git")

## Source functions
sapply(dir("Funktioner",full.names=TRUE), source)

pars = c(theta=10, mu=1, sigma_x=1, sigma_y=1e-2)
N.sim <- 1000
dt.sim <-1e-3
dt.obs <- 1e-2
x0 <- 3
yall <- matrix(data=NA,nrow=1001,ncol=1000)
sda <- matrix(data=NA,nrow = 1000,ncol=1)
for(i in 1:1000){

  
  l <- sim_OU_EM(1000, dt.sim, dt.obs, pars,x0)
    .data <- l$.data
    yall[,i] <- l$x
    sda[i] <- sd(l$x)
}
d_mu <- mean(yall[,1])
sdf <- sd(yall[,1])
d_cov <- acf(yall[,1],pl=F)[1]$acf
d_sig <- 2*sdf*log(sdf/d_cov)
d_theta <- log(sdf/d_cov)

