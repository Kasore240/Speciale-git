library(ggplot2)
library("plotrix")

load("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git/Init_mu_vari.RData")
### plotting for the 20 sim ----

cbind(theta_n,mu_n,sig_x_n,sig_y_n)
theta_n <- a_org[,1:2]
mu_n <- a_org[,3:4]
sig_x_n <- a_org[,5:6]
sig_y_n <- a_org[,7:8]

c <- matrix(data="blue",nrow=21,ncol=1) 
c[21] <- "red"
k <- 1
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



