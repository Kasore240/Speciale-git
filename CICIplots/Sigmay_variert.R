library("plotrix")
setwd("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git")
####---------Make function to Simulate data using Euler Maruyama sim_OU_EM()-----
## Source functions
sapply(dir("Funktioner",full.names=TRUE), source)

imu <- c(2,3)
itheta <- c(1,4)
isx <- c(3,1)
isy <- c(4,2)
int_sig_x <- seq(1e-6,2.5,0.25)
load("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git/data1001/sigmay_varieret_100.RData")
poOLD <- parms_org
sdOLD <- sd_org
load("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git/data1001/TMB_sigmay_varieret_100.RData")
poTMB <- parms_org
sdTMB <- sd_org



#looking at mu
par(mfrow=c(2,1))
CICI(poOLD,imu[2],sdOLD,int_sig_x,rep(1,10))
title("For Old pack",cex.main=0.8)
CICI(poTMB,imu[1],sdTMB,int_sig_x,rep(1,10))
title("For TMB",cex.main=0.8)
mtext("MU when sigma y is changed",side=3,outer=TRUE, line=-1,cex.main=1.5)

#looking at theta
par(mfrow=c(2,1))
CICI(poOLD,itheta[2],sdOLD,int_sig_x,rep(10,10))
title("For Old pack",cex.main=0.8)
CICI(poTMB,itheta[1],sdTMB,int_sig_x,rep(log(10),10))
title("For TMB",cex.main=0.8)
mtext("Theta when sigma y is changed",side=3,outer=TRUE, line=-1,cex.main=1.5)

#lookign at sig x
par(mfrow=c(2,1))
CICI(poOLD,isx[2],sdOLD,int_sig_x,rep(log(1),10))
title("For Old pack",cex.main=0.8)
CICI(poTMB,isx[1],sdTMB,int_sig_x,rep(log(1),10))
title("For TMB",cex.main=0.8)
mtext("Sigma x when sigma y is changed",side=3,outer=TRUE, line=-1,cex.main=1.5)

#looking at sigy
par(mfrow=c(2,1))
CICI(poOLD,isy[2],sdOLD,int_sig_x,log(int_sig_x))
title("For Old pack",cex.main=0.8)
CICI(poTMB,isy[1],sdTMB,int_sig_x,log(int_sig_x))
title("For TMB",cex.main=0.8)
mtext("Sigma y when sigma y is changed",side=3,outer=TRUE, line=-1,cex.main=1.5)


### count nan/na ----
#det er den samme for alle parametre jo.
oi <- 1
ti <- 3
nanssd <- matrix(data = NA, nrow = 10,ncol = 2)
nansp <- nans <- matrix(data = NA, nrow = 10,ncol = 2)

for (i in 1:10){
  nansp[i,1]<- sum(is.na(poOLD[,oi,i]))
  nansp[i,2]<- sum(is.na(poOLD[,ti,i]))
  
  nanssd[i,1] <- sum(is.na(sdOLD[,oi,i]))
  nanssd[i,2] <- sum(is.na(sdTMB[,ti,i]))
  print(100 - sum(is.na(sdTMB[,oi,i])))
  
}
nanssd
nansp
