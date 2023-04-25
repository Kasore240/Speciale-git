library(ggplot2)
library("plotrix")
CICI <- function(po,poi,sd,sdi,xscale,p,expbol){
  confs <- matrix(data = NA, nrow = 10, ncol = 3)
  for (k in 1:10) {
    
    if (expbol) {
      low =exp(c(po[,poi,k]) - (2*c(sd[,sdi,k])))
      up =exp(c(po[,poi,k]) + (2*c(sd[,sdi,k])))
    } else {
      low =(c(po[,poi,k]) - (2*c(sd[,sdi,k])))
      up =(c(po[,poi,k]) + (2*c(sd[,sdi,k])))
    }
    
    x <- sum(p< up & p> low)
    b <- binom.test(x,100,0.95,conf.level = 0.95)
    confs[k,1] <- x/100
    confs[k,2:3] <- b$conf.int
  }
  
  plotCI(x=xscale, y= confs[,1],li=confs[,2],ui=confs[,3])
  abline(h = 0.95, col = "red")
  
}

### plotting for the 20 sim ----

theta_n <- a_org[,1:2]
mu_n <- a_org[,3:4]
sig_x_n <- a_org[,5:6]
sig_y_n <- a_org[,7:8]

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






### confidence of confidence interval MU CTSMR -----
load("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git/data100/Init_mu_vari_100.RData")
confs <- matrix(data = NA, nrow = 10, ncol = 3)
p <- 1
for (k in 1:10) {


low =(c(parms_org[,2,k]) - (2*c(sd_org[,3,k])))
up =(c(parms_org[,2,k]) + (2*c(sd_org[,3,k])))


x <- sum(p< up & p> low)
if (is.na(x)){ x<-100}
b <- binom.test(x,100,0.95,conf.level = 0.95)
confs[k,1] <- x/100
confs[k,2:3] <- b$conf.int

}
int_mu <-  seq(0,4.5,0.5)
plotCI(x=int_mu, y= confs[,1],li=confs[,2],ui=confs[,3])
abline(h = 0.95, col = "red")
line(x=1:10,y=rep(0.95,10))

### confidence of confidence interval MU CTSMRTMB -----
load("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git/data100/TMB_init_mu_varieret_100.RData")
confs <- matrix(data = NA, nrow = 10, ncol = 3)
p <- 1
for (k in 1:10) {
  
  
  low =(c(parms_org[,2,k]) - (2*c(sd_org[,3,k])))
  up =(c(parms_org[,2,k]) + (2*c(sd_org[,3,k])))
  
  
  x <- sum(p< up & p> low)
  b <- binom.test(x,100,0.95,conf.level = 0.95)
  confs[k,1] <- x/100
  confs[k,2:3] <- b$conf.int
  
}
int_mu <-  seq(0,4.5,0.5)
plotCI(x=int_mu, y= confs[,1],li=confs[,2],ui=confs[,3])
abline(h = 0.95, col = "red")

### use new plotting function -----
imu <- c(2,3)
itheta <- c(1,4)
isx <- c(3,1)
isy <- c(4,2)
int_mu <-  seq(0,4.5,0.5)

load("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git/data100/Init_mu_vari_100.RData")
poOLD <- parms_org
sdOLD <- sd_org

load("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git/data100/TMB_init_mu_varieret_100.RData")
poTMB <- parms_org
sdTMB <- sd_org

#looking at mu
par(mfrow=c(2,1))
CICI(poOLD,imu[1],sdOLD,imu[2],int_mu,1,FALSE)
title("For Old pack",cex.main=0.8)
CICI(poTMB,imu[1],sdTMB,imu[2],int_mu,1,FALSE)
title("For TMB",cex.main=0.8)
mtext("MU when initial guess of mu is changed",side=3,outer=TRUE, line=-1,cex.main=1.5)

#looking at theta
par(mfrow=c(2,1))
CICI(poOLD,itheta[1],sdOLD,itheta[2],int_mu,10,FALSE)
title("For Old pack",cex.main=0.8)
CICI(poTMB,itheta[1],sdTMB,itheta[2],int_mu,10,TRUE)
title("For TMB",cex.main=0.8)
mtext("Theta when initial guess of mu is changed",side=3,outer=TRUE, line=-1,cex.main=1.5)

#lookign at sig x
par(mfrow=c(2,1))
CICI(poOLD,isx[1],sdOLD,isx[2],int_mu,1,TRUE)
title("For Old pack",cex.main=0.8)
CICI(poTMB,isx[1],sdTMB,isx[2],int_mu,1,TRUE)
title("For TMB",cex.main=0.8)
mtext("Sigma x when initial guess of mu is changed",side=3,outer=TRUE, line=-1,cex.main=1.5)

#looking at sigy
par(mfrow=c(2,1))
CICI(poOLD,isy[1],sdOLD,isy[2],int_mu,1e-2,TRUE)
title("For Old pack",cex.main=0.8)
CICI(poTMB,isy[1],sdTMB,isy[2],int_mu,1e-2,TRUE)
title("For TMB",cex.main=0.8)
mtext("Sigma y when initial guess of mu is changed",side=3,outer=TRUE, line=-1,cex.main=1.5)

