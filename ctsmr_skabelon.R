# ctsmr ------

#load librarties
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


# call functions ---- 
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

l <- sim_OU_EM(N.sim, dt.sim, dt.obs, pars, x0)
.data <- l$.data
x <- l$x
model <- OU_ctsmr(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub)

#Predict --- 

# Run the parameter estimation
fit <- model$estimate(.data)

# See the summary of the estimation
summ_org <- summary(fit)
summ_org_exp <- c(fit$xm[5], fit$xm[4], exp(fit$xm[2]),(exp(fit$xm[3])))





# ---------------------test influence of sigma_x  ------

theta_n <- matrix(data = NA, nrow=10, ncol =2)
mu_n <- matrix(data = NA, nrow=10, ncol =2)
sig_x_n <- matrix(data = NA, nrow=10, ncol =2)
sig_y_n <- matrix(data = NA, nrow=10, ncol =2)
noise <- seq(1e-1,5,0.49)
sd_mean <- matrix(data = NA, nrow = 10, ncol = 4)
set.seed(11)
parms_org <- array(rep(NaN, 20*4*10), dim=c(20, 4, 10))
sd_org <- array(rep(NaN, 20*4*10), dim=c(20, 4, 10))
N.sim<-20

for(j in 1:10){
pars = c(theta=10, mu=1, sigma_x=noise[j], sigma_y=1e-2)
  
#to easy change the start-parameters c(theta, mu, sigma_x, sigma_y,x0)
init_pars <- c(1, 1.5, 1e-1, 1e-1, 10)

for (i in 1:N.sim){

  l <- sim_OU_EM(1000, dt.sim, dt.obs, pars, x0)
  .data <- l$.data
  x <- l$x
  
  model <- OU_ctsmr(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub)
  
  fit <- model$estimate(.data)
  
  summ_org <- summary(fit)
  summ_org_exp <- c(fit$xm[5], fit$xm[4], exp(fit$xm[2]),(exp(fit$xm[3])))
  parms_org[i,,j] <- summ_org_exp
  tmp <- predict(fit)[[1]]
  sd_org[i,,j] <- fit$sd[2:5]
  
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

  ### plotting for the 20 sim ----
c <- matrix(data="blue",nrow=21,ncol=1) 
c[21] <- "red"
k <- 10
plotCI(x = 1:21,               # plotrix plot with confidence intervals
       y = c(parms_org[,1,k],theta_n[k,2]),
       li = c(parms_org[,1,k],theta_n[k,2]) - (2*c(sd_org[,4,k],sd_mean[k,4])) ,
       ui = c(parms_org[,1,k],theta_n[k,2]) + (2*c(sd_org[,4,k],sd_mean[k,4])),col = c)
title(main= "Theta")

#Mu
plotCI(x = 1:21,               # plotrix plot with confidence intervals
       y = c(parms_org[,2,k],mu_n[k,2]),
       li = c(parms_org[,2,k],mu_n[k,2]) - (2*c(sd_org[,3,k],sd_mean[k,3])) ,
       ui = c(parms_org[,2,k],mu_n[k,2]) + (2*c(sd_org[,3,k],sd_mean[k,3])),col = c)
title(main= "Mu")
#sigma_x 
plotCI(x = 1:21,               # plotrix plot with confidence intervals
       y = c(parms_org[,3,k],sig_x_n[k,2]),
       li = c(parms_org[,3,k],sig_x_n[k,2]) - exp(2*c(sd_org[,1,k],sd_mean[k,1])) ,
       ui = c(parms_org[,3,k],sig_x_n[k,2]) + exp(2*c(sd_org[,1,k],sd_mean[k,1])),col = c)
title(main= "Sigma_x")

#sigma_y
plotCI(x = 1:21,               # plotrix plot with confidence intervals
       y = c(parms_org[,4,k],sig_y_n[k,2]),
       li = c(parms_org[,4,k],sig_y_n[k,2]) - exp(2*c(sd_org[,2,k],sd_mean[k,2])) ,
       ui = c(parms_org[,4,k],sig_y_n[k,2]) + exp(2*c(sd_org[,2,k],sd_mean[k,2])),col = c)
title(main= "Sigma_Y")

    ### plotting after taking mean over the 20 sim ------



#theta
plotCI(x = noise[1:6],               # plotrix plot with confidence intervals
       y = theta_n[1:6,2],
       li = theta_n[1:6,2] - (2*sd_mean[1:6,4]) ,
       ui = theta_n[1:6,2] + (2*sd_mean[1:6,4]))
#Mu
plotCI(x = noise[1:6],               # plotrix plot with confidence intervals
       y = mu_n[1:6,2],
       li = mu_n[1:6,2] - (2*sd_mean[1:6,3]) ,
       ui = mu_n[1:6,2] + (2*sd_mean[1:6,3]))
#sigma x
plotCI(x = noise[1:6],               # plotrix plot with confidence intervals
       y = sig_x_n[1:6,2],
       li = sig_x_n[1:6,2] - exp(2*sd_mean[1:6,1]) ,
       ui = sig_x_n[1:6,2] + exp(2*sd_mean[1:6,1]))
#sigma y
plotCI(x = noise[1:6],               # plotrix plot with confidence intervals
       y = sig_y_n[1:6,2],
       li = sig_y_n[1:6,2] - exp(2*sd_mean[1:6,2]) ,
       ui = sig_y_n[1:6,2] + exp(2*sd_mean[1:6,2]))


# ---------------------test influence of sigma_y (skal lige køres igen med rigtig sig_x) ------

theta_n <- matrix(data = NA, nrow=10, ncol =2)
mu_n <- matrix(data = NA, nrow=10, ncol =2)
sig_x_n <- matrix(data = NA, nrow=10, ncol =2)
sig_y_n <- matrix(data = NA, nrow=10, ncol =2)
noise <-  seq(1e-3,4.5,0.49)
sd_mean <- matrix(data = NA, nrow = 10, ncol = 4)
set.seed(11)
parms_org <- array(rep(NaN, 20*4*10), dim=c(20, 4, 10))
sd_org <- array(rep(NaN, 20*4*10), dim=c(20, 4, 10))
N.sim<-20

for(j in 1:10){
  pars = c(theta=10, mu=1, sigma_x=1, sigma_y=noise[j])
  
  #to easy change the start-parameters c(theta, mu, sigma_x, sigma_y,x0)
  init_pars <- c(1, 1.5, 1e-1, 1e-1, 10)
  
  for (i in 1:N.sim){
    
    l <- sim_OU_EM(1000, dt.sim, dt.obs, pars, x0)
    .data <- l$.data
    x <- l$x
    
    model <- OU_ctsmr(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub)
    
    fit <- model$estimate(.data)
    
    summ_org <- summary(fit)
    summ_org_exp <- c(fit$xm[5], fit$xm[4], exp(fit$xm[2]),(exp(fit$xm[3])))
    parms_org[i,,j] <- summ_org_exp
    tmp <- predict(fit)[[1]]
    sd_org[i,,j] <- fit$sd[2:5]
    
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
save(a_org,parms_org,sd_org,sd_mean,file = "sigmay_varieret.RData")

  ### plotting for the 20 sim ----
c <- matrix(data="blue",nrow=21,ncol=1) 
c[21] <- "red"
  k <- 5
  plotCI(x = 1:21,               # plotrix plot with confidence intervals
         y = c(parms_org[,1,k],theta_n[k,2]),
         li = c(parms_org[,1,k],theta_n[k,2]) - (2*c(sd_org[,4,k],sd_mean[k,4])) ,
         ui = c(parms_org[,1,k],theta_n[k,2]) + (2*c(sd_org[,4,k],sd_mean[k,4])),col = c)
  title(main= "Theta")
  
  #Mu
  plotCI(x = 1:21,               # plotrix plot with confidence intervals
         y = c(parms_org[,2,k],mu_n[k,2]),
         li = c(parms_org[,2,k],mu_n[k,2]) - (2*c(sd_org[,3,k],sd_mean[k,3])) ,
         ui = c(parms_org[,2,k],mu_n[k,2]) + (2*c(sd_org[,3,k],sd_mean[k,3])),col = c)
  title(main= "Mu")
  #sigma_x 
  plotCI(x = 1:21,               # plotrix plot with confidence intervals
         y = c(parms_org[,3,k],sig_x_n[k,2]),
         li = c(parms_org[,3,k],sig_x_n[k,2]) - exp(2*c(sd_org[,1,k],sd_mean[k,1])) ,
         ui = c(parms_org[,3,k],sig_x_n[k,2]) + exp(2*c(sd_org[,1,k],sd_mean[k,1])),col = c)
  title(main= "Sigma_x")
  
  #sigma_y
  plotCI(x = 1:21,               # plotrix plot with confidence intervals
         y = c(parms_org[,4,k],sig_y_n[k,2]),
         li = c(parms_org[,4,k],sig_y_n[k,2]) - (2*c(sd_org[,2,k],sd_mean[k,2])) ,
         ui = c(parms_org[,4,k],sig_y_n[k,2]) + (2*c(sd_org[,2,k],sd_mean[k,2])),col = c)
  title(main= "Sigma_Y")
  
    ### plotting after taking mean over the 20 sim ------
  
  
  
  #theta
  plotCI(x = noise[1:6],               # plotrix plot with confidence intervals
         y = theta_n[1:6,2],
         li = theta_n[1:6,2] - (2*sd_mean[1:6,4]) ,
         ui = theta_n[1:6,2] + (2*sd_mean[1:6,4]))
  #Mu
  plotCI(x = noise[1:6],               # plotrix plot with confidence intervals
         y = mu_n[1:6,2],
         li = mu_n[1:6,2] - (2*sd_mean[1:6,3]) ,
         ui = mu_n[1:6,2] + (2*sd_mean[1:6,3]))
  #sigma x
  plotCI(x = noise[1:6],               # plotrix plot with confidence intervals
         y = sig_x_n[1:6,2],
         li = sig_x_n[1:6,2] - exp(2*sd_mean[1:6,1]) ,
         ui = sig_x_n[1:6,2] + exp(2*sd_mean[1:6,1]))
  #sigma y
  plotCI(x = noise[1:6],               # plotrix plot with confidence intervals
         y = sig_y_n[1:6,2],
         li = sig_y_n[1:6,2] - exp(2*sd_mean[1:6,2]) ,
         ui = sig_y_n[1:6,2] + exp(2*sd_mean[1:6,2]))
  
  
  
# --------------------test influence of mu start guess  ------
  
  theta_n <- matrix(data = NA, nrow=10, ncol =2)
  mu_n <- matrix(data = NA, nrow=10, ncol =2)
  sig_x_n <- matrix(data = NA, nrow=10, ncol =2)
  sig_y_n <- matrix(data = NA, nrow=10, ncol =2)
  int_mu <-  seq(0,4.5,0.5)
  sd_mean <- matrix(data = NA, nrow = 10, ncol = 4)
  set.seed(11)
  parms_org <- array(rep(NaN, 20*4*10), dim=c(20, 4, 10))
  sd_org <- array(rep(NaN, 20*4*10), dim=c(20, 4, 10))
  N.sim<-20
  
  #to easy change lower bounds  c(theta, mu, sigma_x, sigma_y,x0)
  init_lb <- c(-1e-5, -1e-5, 1e-10, 1e-12, -1)
  
  #to easy change upper bounds  c(theta, mu, sigma_x, sigma_y,x0)
  init_ub <- c(50, 5, 10, 10, 100)
  
  for(j in 1:10){
    pars = c(theta=10, mu=1, sigma_x=1, sigma_y=1e-2)
    
    #to easy change the start-parameters c(theta, mu, sigma_x, sigma_y,x0)
    init_pars <- c(1, int_mu[j], 1e-1, 1e-1, 10)
    
    for (i in 1:N.sim){
      
      l <- sim_OU_EM(1000, dt.sim, dt.obs, pars, init_pars[5])
      .data <- l$.data
      x <- l$x
      
      model <- OU_ctsmr(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub)
      
      fit <- model$estimate(.data)
      
      summ_org <- summary(fit)
      summ_org_exp <- c(fit$xm[5], fit$xm[4], fit$xm[2],fit$xm[3])
      parms_org[i,,j] <- summ_org_exp
      tmp <- predict(fit)[[1]]
      sd_org[i,,j] <- fit$sd[2:5]
      
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
  
  ### plotting for the 20 sim ----
  c <- matrix(data="blue",nrow=21,ncol=1) 
  c[21] <- "red"
    k <- 5
    plotCI(x = 1:21,               # plotrix plot with confidence intervals
           y = c(parms_org[,1,k],theta_n[k,2]),
           li = c(parms_org[,1,k],theta_n[k,2]) - (2*c(sd_org[,4,k],sd_mean[k,4])) ,
           ui = c(parms_org[,1,k],theta_n[k,2]) + (2*c(sd_org[,4,k],sd_mean[k,4])),col = c)
    title(main= "Theta")
    
    #Mu
    plotCI(x = 1:21,               # plotrix plot with confidence intervals
           y = c(parms_org[,2,k],mu_n[k,2]),
           li = c(parms_org[,2,k],mu_n[k,2]) - (2*c(sd_org[,3,k],sd_mean[k,3])) ,
           ui = c(parms_org[,2,k],mu_n[k,2]) + (2*c(sd_org[,3,k],sd_mean[k,3])),col = c)
    title(main= "Mu")
    #sigma_x 
    plotCI(x = 1:21,               # plotrix plot with confidence intervals
           y = exp(c(parms_org[,3,k],sig_x_n[k,2])),
           li = exp(c(parms_org[,3,k],sig_x_n[k,2]) - 2*c(sd_org[,1,k],sd_mean[k,1])) ,
           ui = exp(c(parms_org[,3,k],sig_x_n[k,2]) + 2*c(sd_org[,1,k],sd_mean[k,1])),col = c)
    title(main= "Sigma_x")
    
    #sigma_y
    plotCI(x = 1:21,               # plotrix plot with confidence intervals
           y = (c(parms_org[,4,k],sig_y_n[k,2])),
           li = (c(parms_org[,4,k],sig_y_n[k,2]) - 2*c(sd_org[,2,k],sd_mean[k,2])) ,
           ui = (c(parms_org[,4,k],sig_y_n[k,2]) + 2*c(sd_org[,2,k],sd_mean[k,2])),col = c)
    title(main= "Sigma_Y")
    
    ### plotting after taking mean over the 20 sim ------
    
    
    
    #theta
    plotCI(x = noise[1:6],               # plotrix plot with confidence intervals
           y = theta_n[1:6,2],
           li = theta_n[1:6,2] - (2*sd_mean[1:6,4]) ,
           ui = theta_n[1:6,2] + (2*sd_mean[1:6,4]))
    #Mu
    plotCI(x = noise[1:6],               # plotrix plot with confidence intervals
           y = mu_n[1:6,2],
           li = mu_n[1:6,2] - (2*sd_mean[1:6,3]) ,
           ui = mu_n[1:6,2] + (2*sd_mean[1:6,3]))
    #sigma x
    plotCI(x = noise[1:6],               # plotrix plot with confidence intervals
           y = sig_x_n[1:6,2],
           li = exp(sig_x_n[1:6,2] - 2*sd_mean[1:6,1]) ,
           ui = exp(sig_x_n[1:6,2] + 2*sd_mean[1:6,1]))
    #sigma y
    plotCI(x = noise[1:6],               # plotrix plot with confidence intervals
           y = sig_y_n[1:6,2],
           li = sig_y_n[1:6,2] - exp(2*sd_mean[1:6,2]) ,
           ui = sig_y_n[1:6,2] + exp(2*sd_mean[1:6,2]))
    
  
    
    

# --------------------test influence of theta start guess  ------
    
    theta_n <- matrix(data = NA, nrow=10, ncol =2)
    mu_n <- matrix(data = NA, nrow=10, ncol =2)
    sig_x_n <- matrix(data = NA, nrow=10, ncol =2)
    sig_y_n <- matrix(data = NA, nrow=10, ncol =2)
    int_theta <-  seq(0,18,2)
    sd_mean <- matrix(data = NA, nrow = 10, ncol = 4)
    set.seed(11)
    parms_org <- array(rep(NaN, 20*4*10), dim=c(20, 4, 10))
    sd_org <- array(rep(NaN, 20*4*10), dim=c(20, 4, 10))
    N.sim<-20
    
    #to easy change lower bounds  c(theta, mu, sigma_x, sigma_y,x0)
    init_lb <- c(-1e-5, 0, 1e-10, 1e-12, -1)
    
    #to easy change upper bounds  c(theta, mu, sigma_x, sigma_y,x0)
    init_ub <- c(50, 5, 10, 10, 100)
    
    for(j in 1:10){
      pars = c(theta=10, mu=1, sigma_x=1, sigma_y=1e-2)
      
      #to easy change the start-parameters c(theta, mu, sigma_x, sigma_y,x0)
      init_pars <- c(int_theta[j], 1.5, 1e-1, 1e-1, 10)
      
      for (i in 1:N.sim){
        
        l <- sim_OU_EM(1000, dt.sim, dt.obs, pars, x0)
        .data <- l$.data
        x <- l$x
        
        model <- OU_ctsmr(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub)
        
        fit <- model$estimate(.data)
        
        summ_org <- summary(fit)
        summ_org_exp <- c(fit$xm[5], fit$xm[4], exp(fit$xm[2]),(exp(fit$xm[3])))
        parms_org[i,,j] <- summ_org_exp
        tmp <- predict(fit)[[1]]
        sd_org[i,,j] <- fit$sd[2:5]
        
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
    
  ### plotting for the 20 sim ----
    c <- matrix(data="blue",nrow=21,ncol=1) 
    c[21] <- "red"
      k <- 10
      plotCI(x = 1:21,               # plotrix plot with confidence intervals
             y = c(parms_org[,1,k],theta_n[k,2]),
             li = c(parms_org[,1,k],theta_n[k,2]) - (2*c(sd_org[,4,k],sd_mean[k,4])) ,
             ui = c(parms_org[,1,k],theta_n[k,2]) + (2*c(sd_org[,4,k],sd_mean[k,4])),col = c)
      title(main= "Theta")
      
      #Mu
      plotCI(x = 1:21,               # plotrix plot with confidence intervals
             y = c(parms_org[,2,k],mu_n[k,2]),
             li = c(parms_org[,2,k],mu_n[k,2]) - (2*c(sd_org[,3,k],sd_mean[k,3])) ,
             ui = c(parms_org[,2,k],mu_n[k,2]) + (2*c(sd_org[,3,k],sd_mean[k,3])),col = c)
      title(main= "Mu")
      #sigma_x 
      plotCI(x = 1:21,               # plotrix plot with confidence intervals
             y = c(parms_org[,3,k],sig_x_n[k,2]),
             li = c(parms_org[,3,k],sig_x_n[k,2]) - exp(2*c(sd_org[,1,k],sd_mean[k,1])) ,
             ui = c(parms_org[,3,k],sig_x_n[k,2]) + exp(2*c(sd_org[,1,k],sd_mean[k,1])),col = c)
      title(main= "Sigma_x")
      
      #sigma_y
      plotCI(x = 1:21,               # plotrix plot with confidence intervals
             y = c(parms_org[,4,k],sig_y_n[k,2]),
             li = c(parms_org[,4,k],sig_y_n[k,2]) - (2*c(sd_org[,2,k],sd_mean[k,2])) ,
             ui = c(parms_org[,4,k],sig_y_n[k,2]) + (2*c(sd_org[,2,k],sd_mean[k,2])),col = c)
      title(main= "Sigma_Y")
      
    ### plotting after taking mean over the 20 sim ------
      
      
      
      #theta
      plotCI(x = int_theta,               # plotrix plot with confidence intervals
             y = theta_n[,2],
             li = theta_n[,2] - (2*sd_mean[,4]) ,
             ui = theta_n[,2] + (2*sd_mean[,4]))
      #Mu
      plotCI(x = int_theta,               # plotrix plot with confidence intervals
             y = mu_n[,2],
             li = mu_n[,2] - (2*sd_mean[,3]) ,
             ui = mu_n[,2] + (2*sd_mean[,3]))
      #sigma x
      plotCI(x = int_theta,               # plotrix plot with confidence intervals
             y = sig_x_n[,2],
             li = sig_x_n[,2] - exp(2*sd_mean[,1]) ,
             ui = sig_x_n[,2] + exp(2*sd_mean[,1]))
      #sigma y
      plotCI(x = int_theta,               # plotrix plot with confidence intervals
             y = sig_y_n[,2],
             li = sig_y_n[,2] - (2*sd_mean[,2]) ,
             ui = sig_y_n[,2] + (2*sd_mean[,2]))
      
      
    
    
    

      
#--------------------test influence of x0 start guess  (ikke kørt og gemt )------
      
      theta_n <- matrix(data = NA, nrow=10, ncol =2)
      mu_n <- matrix(data = NA, nrow=10, ncol =2)
      sig_x_n <- matrix(data = NA, nrow=10, ncol =2)
      sig_y_n <- matrix(data = NA, nrow=10, ncol =2)
      int_x0 <-  seq(-1,12.5,1.5)
      sd_mean <- matrix(data = NA, nrow = 10, ncol = 4)
      set.seed(11)
      parms_org <- array(rep(NaN, 20*4*10), dim=c(20, 4, 10))
      sd_org <- array(rep(NaN, 20*4*10), dim=c(20, 4, 10))
      N.sim<-20
      
      #to easy change lower bounds  c(theta, mu, sigma_x, sigma_y,x0)
      init_lb <- c(-1e-5, 0, 1e-10, 1e-12, -1)
      
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
          summ_org_exp <- c(fit$xm[5], fit$xm[4], exp(fit$xm[2]),(exp(fit$xm[3])))
          parms_org[i,,j] <- summ_org_exp
          tmp <- predict(fit)[[1]]
          sd_org[i,,j] <- fit$sd[2:5]
          
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
      save(a_org,parms_org,sd_org,sd_mean,file = "init_x0_varieret.RData")
      
  ### plotting for the 20 sim ----
      c <- matrix(data="blue",nrow=21,ncol=1) 
      c[21] <- "red"
        k <- 10
        plotCI(x = 1:21,               # plotrix plot with confidence intervals
               y = c(parms_org[,1,k],theta_n[k,2]),
               li = c(parms_org[,1,k],theta_n[k,2]) - (2*c(sd_org[,4,k],sd_mean[k,4])) ,
               ui = c(parms_org[,1,k],theta_n[k,2]) + (2*c(sd_org[,4,k],sd_mean[k,4])),col = c)
        title(main= "Theta")
        
        #Mu
        plotCI(x = 1:21,               # plotrix plot with confidence intervals
               y = c(parms_org[,2,k],mu_n[k,2]),
               li = c(parms_org[,2,k],mu_n[k,2]) - (2*c(sd_org[,3,k],sd_mean[k,3])) ,
               ui = c(parms_org[,2,k],mu_n[k,2]) + (2*c(sd_org[,3,k],sd_mean[k,3])),col = c)
        title(main= "Mu")
        #sigma_x 
        plotCI(x = 1:21,               # plotrix plot with confidence intervals
               y = c(parms_org[,3,k],sig_x_n[k,2]),
               li = c(parms_org[,3,k],sig_x_n[k,2]) - exp(2*c(sd_org[,1,k],sd_mean[k,1])) ,
               ui = c(parms_org[,3,k],sig_x_n[k,2]) + exp(2*c(sd_org[,1,k],sd_mean[k,1])),col = c)
        title(main= "Sigma_x")
        
        #sigma_y
        plotCI(x = 1:21,               # plotrix plot with confidence intervals
               y = c(parms_org[,4,k],sig_y_n[k,2]),
               li = c(parms_org[,4,k],sig_y_n[k,2]) - (2*c(sd_org[,2,k],sd_mean[k,2])) ,
               ui = c(parms_org[,4,k],sig_y_n[k,2]) + (2*c(sd_org[,2,k],sd_mean[k,2])),col = c)
        title(main= "Sigma_Y")
        
    ### plotting after taking mean over the 20 sim ------
        
        
        
        #theta
        plotCI(x = int_theta,               # plotrix plot with confidence intervals
               y = theta_n[,2],
               li = theta_n[,2] - (2*sd_mean[,4]) ,
               ui = theta_n[,2] + (2*sd_mean[,4]))
        #Mu
        plotCI(x = int_theta,               # plotrix plot with confidence intervals
               y = mu_n[,2],
               li = mu_n[,2] - (2*sd_mean[,3]) ,
               ui = mu_n[,2] + (2*sd_mean[,3]))
        #sigma x
        plotCI(x = int_theta,               # plotrix plot with confidence intervals
               y = sig_x_n[,2],
               li = sig_x_n[,2] - exp(2*sd_mean[,1]) ,
               ui = sig_x_n[,2] + exp(2*sd_mean[,1]))
        #sigma y
        plotCI(x = int_theta,               # plotrix plot with confidence intervals
               y = sig_y_n[,2],
               li = sig_y_n[,2] - (2*sd_mean[,2]) ,
               ui = sig_y_n[,2] + (2*sd_mean[,2]))
        
        
        
        
# -------------------test influence of sig_x start guess  (ikke kørt og gemt )------
        
theta_n <- matrix(data = NA, nrow=10, ncol =2)
mu_n <- matrix(data = NA, nrow=10, ncol =2)
sig_x_n <- matrix(data = NA, nrow=10, ncol =2)
sig_y_n <- matrix(data = NA, nrow=10, ncol =2)
int_sig_x <-  seq(1e-10,3,0.3)
sd_mean <- matrix(data = NA, nrow = 10, ncol = 4)
set.seed(11)
parms_org <- array(rep(NaN, 20*4*10), dim=c(20, 4, 10))
sd_org <- array(rep(NaN, 20*4*10), dim=c(20, 4, 10))
N.sim<-20

#to easy change lower bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_lb <- c(-1e-5, 0, 1e-10, 1e-12, -1)

#to easy change upper bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_ub <- c(50, 5, 10, 10, 100)

for(j in 1:10){
  pars = c(theta=10, mu=1, sigma_x=1, sigma_y=1e-2)
  
  #to easy change the start-parameters c(theta, mu, sigma_x, sigma_y,x0)
  init_pars <- c(1, 1.5, int_sig_x[j], 1e-1, 10)
  
  for (i in 1:N.sim){
    
    l <- sim_OU_EM(1000, dt.sim, dt.obs, pars, 3)
    .data <- l$.data
    x <- l$x
    
    model <- OU_ctsmr(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub)
    
    fit <- model$estimate(.data)
    
    summ_org <- summary(fit)
    summ_org_exp <- c(fit$xm[5], fit$xm[4], exp(fit$xm[2]),(exp(fit$xm[3])))
    parms_org[i,,j] <- summ_org_exp
    tmp <- predict(fit)[[1]]
    sd_org[i,,j] <- fit$sd[2:5]
    
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

save(a_org,parms_org,sd_org,sd_mean,file = "init_sigmax_varieret.RData")

  ### plotting for the 20 sim ----
        c <- matrix(data="blue",nrow=21,ncol=1) 
        c[21] <- "red"
          k <- 10
          plotCI(x = 1:21,               # plotrix plot with confidence intervals
                 y = c(parms_org[,1,k],theta_n[k,2]),
                 li = c(parms_org[,1,k],theta_n[k,2]) - (2*c(sd_org[,4,k],sd_mean[k,4])) ,
                 ui = c(parms_org[,1,k],theta_n[k,2]) + (2*c(sd_org[,4,k],sd_mean[k,4])),col = c)
          title(main= "Theta")
          
          #Mu
          plotCI(x = 1:21,               # plotrix plot with confidence intervals
                 y = c(parms_org[,2,k],mu_n[k,2]),
                 li = c(parms_org[,2,k],mu_n[k,2]) - (2*c(sd_org[,3,k],sd_mean[k,3])) ,
                 ui = c(parms_org[,2,k],mu_n[k,2]) + (2*c(sd_org[,3,k],sd_mean[k,3])),col = c)
          title(main= "Mu")
          #sigma_x 
          plotCI(x = 1:21,               # plotrix plot with confidence intervals
                 y = c(parms_org[,3,k],sig_x_n[k,2]),
                 li = c(parms_org[,3,k],sig_x_n[k,2]) - exp(2*c(sd_org[,1,k],sd_mean[k,1])) ,
                 ui = c(parms_org[,3,k],sig_x_n[k,2]) + exp(2*c(sd_org[,1,k],sd_mean[k,1])),col = c)
          title(main= "Sigma_x")
          
          #sigma_y
          plotCI(x = 1:21,               # plotrix plot with confidence intervals
                 y = c(parms_org[,4,k],sig_y_n[k,2]),
                 li = c(parms_org[,4,k],sig_y_n[k,2]) - (2*c(sd_org[,2,k],sd_mean[k,2])) ,
                 ui = c(parms_org[,4,k],sig_y_n[k,2]) + (2*c(sd_org[,2,k],sd_mean[k,2])),col = c)
          title(main= "Sigma_Y")
          
    ### plotting after taking mean over the 20 sim ------
          
          
          
          #theta
          plotCI(x = int_theta,               # plotrix plot with confidence intervals
                 y = theta_n[,2],
                 li = theta_n[,2] - (2*sd_mean[,4]) ,
                 ui = theta_n[,2] + (2*sd_mean[,4]))
          #Mu
          plotCI(x = int_theta,               # plotrix plot with confidence intervals
                 y = mu_n[,2],
                 li = mu_n[,2] - (2*sd_mean[,3]) ,
                 ui = mu_n[,2] + (2*sd_mean[,3]))
          #sigma x
          plotCI(x = int_theta,               # plotrix plot with confidence intervals
                 y = sig_x_n[,2],
                 li = sig_x_n[,2] - exp(2*sd_mean[,1]) ,
                 ui = sig_x_n[,2] + exp(2*sd_mean[,1]))
          #sigma y
          plotCI(x = int_theta,               # plotrix plot with confidence intervals
                 y = sig_y_n[,2],
                 li = sig_y_n[,2] - (2*sd_mean[,2]) ,
                 ui = sig_y_n[,2] + (2*sd_mean[,2]))
          
        
        
        
      
      

      
          

# -------------------test influence of sig_y start guess  (ikke kørt og gemt )------

theta_n <- matrix(data = NA, nrow=10, ncol =2)
mu_n <- matrix(data = NA, nrow=10, ncol =2)
sig_x_n <- matrix(data = NA, nrow=10, ncol =2)
sig_y_n <- matrix(data = NA, nrow=10, ncol =2)
int_sig_y <-  seq(1e-5,1,1e-1)
sd_mean <- matrix(data = NA, nrow = 10, ncol = 4)
set.seed(11)
parms_org <- array(rep(NaN, 20*4*10), dim=c(20, 4, 10))
sd_org <- array(rep(NaN, 20*4*10), dim=c(20, 4, 10))
N.sim<-20

#to easy change lower bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_lb <- c(-1e-5, 0, 1e-10, 1e-12, -1)

#to easy change upper bounds  c(theta, mu, sigma_x, sigma_y,x0)
init_ub <- c(50, 5, 10, 10, 100)

for(j in 1:10){
  pars = c(theta=10, mu=1, sigma_x=1, sigma_y=1e-2)
  
  #to easy change the start-parameters c(theta, mu, sigma_x, sigma_y,x0)
  init_pars <- c(1, 1.5, 1e-1, int_sig_y[j],10)
  
  for (i in 1:N.sim){
    
    l <- sim_OU_EM(1000, dt.sim, dt.obs, pars, 3)
    .data <- l$.data
    x <- l$x
    
    model <- OU_ctsmr(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub)
    
    fit <- model$estimate(.data)
    
    summ_org <- summary(fit)
    summ_org_exp <- c(fit$xm[5], fit$xm[4], exp(fit$xm[2]),(exp(fit$xm[3])))
    parms_org[i,,j] <- summ_org_exp
    tmp <- predict(fit)[[1]]
    sd_org[i,,j] <- fit$sd[2:5]
    
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
save(a_org,parms_org,sd_org,sd_mean,file = "init_sigmay_varieret.RData")
          
  ### plotting for the 20 sim ----
          c <- matrix(data="blue",nrow=21,ncol=1) 
          c[21] <- "red"
            k <- 10
            plotCI(x = 1:21,               # plotrix plot with confidence intervals
                   y = c(parms_org[,1,k],theta_n[k,2]),
                   li = c(parms_org[,1,k],theta_n[k,2]) - (2*c(sd_org[,4,k],sd_mean[k,4])) ,
                   ui = c(parms_org[,1,k],theta_n[k,2]) + (2*c(sd_org[,4,k],sd_mean[k,4])),col = c)
            title(main= "Theta")
            
            #Mu
            plotCI(x = 1:21,               # plotrix plot with confidence intervals
                   y = c(parms_org[,2,k],mu_n[k,2]),
                   li = c(parms_org[,2,k],mu_n[k,2]) - (2*c(sd_org[,3,k],sd_mean[k,3])) ,
                   ui = c(parms_org[,2,k],mu_n[k,2]) + (2*c(sd_org[,3,k],sd_mean[k,3])),col = c)
            title(main= "Mu")
            #sigma_x 
            plotCI(x = 1:21,               # plotrix plot with confidence intervals
                   y = c(parms_org[,3,k],sig_x_n[k,2]),
                   li = c(parms_org[,3,k],sig_x_n[k,2]) - exp(2*c(sd_org[,1,k],sd_mean[k,1])) ,
                   ui = c(parms_org[,3,k],sig_x_n[k,2]) + exp(2*c(sd_org[,1,k],sd_mean[k,1])),col = c)
            title(main= "Sigma_x")
            
            #sigma_y
            plotCI(x = 1:21,               # plotrix plot with confidence intervals
                   y = c(parms_org[,4,k],sig_y_n[k,2]),
                   li = c(parms_org[,4,k],sig_y_n[k,2]) - (2*c(sd_org[,2,k],sd_mean[k,2])) ,
                   ui = c(parms_org[,4,k],sig_y_n[k,2]) + (2*c(sd_org[,2,k],sd_mean[k,2])),col = c)
            title(main= "Sigma_Y")
            
    ### plotting after taking mean over the 20 sim ------
      #theta
      plotCI(x = int_theta,               # plotrix plot with confidence intervals
             y = theta_n[,2],
             li = theta_n[,2] - (2*sd_mean[,4]) ,
             ui = theta_n[,2] + (2*sd_mean[,4]))
      #Mu
      plotCI(x = int_theta,               # plotrix plot with confidence intervals
             y = mu_n[,2],
             li = mu_n[,2] - (2*sd_mean[,3]) ,
             ui = mu_n[,2] + (2*sd_mean[,3]))
      #sigma x
      plotCI(x = int_theta,               # plotrix plot with confidence intervals
             y = sig_x_n[,2],
             li = sig_x_n[,2] - exp(2*sd_mean[,1]) ,
             ui = sig_x_n[,2] + exp(2*sd_mean[,1]))
      #sigma y
      plotCI(x = int_theta,               # plotrix plot with confidence intervals
             y = sig_y_n[,2],
             li = sig_y_n[,2] - (2*sd_mean[,2]) ,
             ui = sig_y_n[,2] + (2*sd_mean[,2]))
            
      
      
      
# random testing testing ----
  pars = c(theta=10, mu=1, sigma_x=1, sigma_y=1e-2)
  
  #to easy change the start-parameters c(theta, mu, sigma_x, sigma_y,x0)
  init_pars <- c(1, 7.9, 1e-1, 1e-1, 10)
  
  #to easy change lower bounds  c(theta, mu, sigma_x, sigma_y,x0)
  init_lb <- c(-1e-5, 0, 1e-10, 1e-12, -1)
  
  #to easy change upper bounds  c(theta, mu, sigma_x, sigma_y,x0)
  init_ub <- c(50, 8, 10, 10, 100)
  
  N.sim <- 1000
  dt.sim <-1e-3
  dt.obs <- 1e-2
  
  l <- sim_OU_EM(N.sim, dt.sim, dt.obs, pars, 10)
  .data <- l$.data
  x <- l$x
  model <- OU_ctsmr(init_pars=init_pars,init_lb=init_lb, init_ub=init_ub)
  
  #Predict --- 
  
  # Run the parameter estimation
  fit <- model$estimate(.data)
  
  # See the summary of the estimation
  summ_org4<- summary(fit, extended = TRUE)
  summ_org_exp4 <- c(fit$xm[5], fit$xm[4], exp(fit$xm[2]),(exp(fit$xm[3])))
  