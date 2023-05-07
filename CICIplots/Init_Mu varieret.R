
library("plotrix")
CICI <- function(po,index,sd,xscale,p){
  confs <- matrix(data = NA, nrow = 10, ncol = 3)
  for (k in 1:10) {
      low =(c(po[,index,k]) - (2*c(sd[,index,k])))
      up =(c(po[,index,k]) + (2*c(sd[,index,k])))
    x <- sum(p< up & p> low,na.rm=T)
    b <- binom.test(x,100,0.95,conf.level = 0.95)
    confs[k,1] <- x/100
    confs[k,2:3] <- b$conf.int
  }
  
  plotCI(x=xscale, y= confs[,1],li=confs[,2],ui=confs[,3])
  abline(h = 0.95, col = "red")
  
}

imu <- c(2,3)
itheta <- c(1,4)
isx <- c(3,1)
isy <- c(4,2)
int_mu <-  seq(0,4.5,0.5)

load("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git/data1001/Init_mu_vari_100.RData")
poOLD <- parms_org
sdOLD <- sd_org

load("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git/data1001/TMB_init_mu_varieret_100.RData")
poTMB <- parms_org
sdTMB <- sd_org

#looking at mu
par(mfrow=c(2,1))
CICI(poOLD,imu[2],sdOLD,int_mu,1)
title("For Old pack",cex.main=0.8)
CICI(poTMB,imu[1],sdTMB,int_mu,1)
title("For TMB",cex.main=0.8)
mtext("MU when initial guess of mu is changed",side=3,outer=TRUE, line=-1,cex.main=1.5)

#looking at theta
par(mfrow=c(2,1))
CICI(poOLD,itheta[2],sdOLD,int_mu,10)
title("For Old pack",cex.main=0.8)
CICI(poTMB,itheta[1],sdTMB,int_mu,log(10))
title("For TMB",cex.main=0.8)
mtext("Theta when initial guess of mu is changed",side=3,outer=TRUE, line=-1,cex.main=1.5)

#lookign at sig x
par(mfrow=c(2,1))
CICI(poOLD,isx[2],sdOLD,int_mu,log(1))
title("For Old pack",cex.main=0.8)
CICI(poTMB,isx[1],sdTMB,int_mu,log(1))
title("For TMB",cex.main=0.8)
mtext("Sigma x when initial guess of mu is changed",side=3,outer=TRUE, line=-1,cex.main=1.5)

#looking at sigy
par(mfrow=c(2,1))
CICI(poOLD,isy[2],sdOLD,int_mu,log(1e-2))
title("For Old pack",cex.main=0.8)
CICI(poTMB,isy[1],sdTMB,int_mu,log(1e-2))
title("For TMB",cex.main=0.8)
mtext("Sigma y when initial guess of mu is changed",side=3,outer=TRUE, line=-1,cex.main=1.5)

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
  
}
nanssd
nansp
