library("plotrix")
CICI <- function(po,poi,sd,sdi,xscale,p,expbol){
  confs <- matrix(data = NA, nrow = 10, ncol = 3)
  for (k in 1:10) {
    
    if (expbol) {
      low =exp(po[,poi,k] - (2*sd[,sdi,k]))
      up =exp(po[,poi,k] + (2*sd[,sdi,k]))
    } else {
      low =po[,poi,k] - (2*sd[,sdi,k])
      up =po[,poi,k] + (2*sd[,sdi,k])
    }
    
    x <- sum(p< up & p> low)
    if (is.na(x)){ x<-100}
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
int_sig_x <-  seq(1e-10,3,0.3)

load("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git/data100/init_sigmax_varieret_100.RData")
poOLD <- parms_org
sdOLD <- sd_org

load("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git/data100/TMB_init_sig_x_varieret_100.RData")
poTMB <- parms_org
sdTMB <- sd_org

#looking at mu
par(mfrow=c(2,1))
CICI(poOLD,imu[1],sdOLD,imu[2],int_sig_x,1,FALSE)
title("For Old pack",cex.main=0.8)
CICI(poTMB,imu[1],sdTMB,imu[1],int_sig_x,1,FALSE)
title("For TMB",cex.main=0.8)
mtext("MU when initial guess of sigma x is changed",side=3,outer=TRUE, line=-1,cex.main=1.5)

#looking at theta
par(mfrow=c(2,1))
CICI(poOLD,itheta[2],sdOLD,itheta[2],int_sig_x,10,FALSE)
title("For Old pack",cex.main=0.8)
CICI(poTMB,itheta[1],sdTMB,itheta[1],int_sig_x,log(10),FALSE)
title("For TMB",cex.main=0.8)
mtext("Theta when initial guess of sigma x is changed",side=3,outer=TRUE, line=-1,cex.main=1.5)

#lookign at sig x
par(mfrow=c(2,1))
CICI(poOLD,isx[2],sdOLD,isx[2],int_sig_x,log(1),FALSE)
title("For Old pack",cex.main=0.8)
CICI(poTMB,isx[1],sdTMB,isx[1],int_sig_x,log(1),FALSE)
title("For TMB",cex.main=0.8)
mtext("Sigma x when initial guess of sigma x is changed",side=3,outer=TRUE, line=-1,cex.main=1.5)

#looking at sigy
par(mfrow=c(2,1))
CICI(poOLD,isy[2],sdOLD,isy[2],int_sig_x,log(1e-2),FALSE)
title("For Old pack",cex.main=0.8)
CICI(poTMB,isy[1],sdTMB,isy[1],int_sig_x,log(1e-2),FALSE)
title("For TMB",cex.main=0.8)
mtext("Sigma y when initial guess of sigma x is changed",side=3,outer=TRUE, line=-1,cex.main=1.5)



#### test ----------

confs <- matrix(data = NA, nrow = 10, ncol = 3)

    low =exp(c(poTMB[,4,k]) - (2*c(sdTMB[,2,k])))
    up =exp(c(poTMB[,4,k]) + (2*c(sdTMB[,2,k])))


  x <- sum(p< up & p> low)
  if (is.na(x)){ x<-100}
  b <- binom.test(x,100,0.95,conf.level = 0.95)
  confs[k,1] <- x/100
  confs[k,2:3] <- b$conf.int

plotCI(x=xscale, y= confs[,1],li=confs[,2],ui=confs[,3])
abline(h = 0.95, col = "red")
