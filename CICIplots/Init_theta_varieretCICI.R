library("plotrix")
CICI <- function(po,index,sd,xscale,p,v){
  confs <- matrix(data = NA, nrow = 10, ncol = 3)
  for (k in 1:10) {
    pol <- po[sd[,index,k]<v,index,k]
    sdl <- sd[sd[,index,k]<v,index,k]
    low =(c(pol) - (2*c(sdl)))
    up =(c(pol) + (2*c(sdl)))
    x <- sum(p< up & p> low,na.rm=T)
    
    b <- binom.test(x,length(low),0.95,conf.level = 0.95)
    confs[k,1] <- x/100
    confs[k,2:3] <- b$conf.int
  }
  
  plotCI(x=xscale, y= confs[,1],li=confs[,2],ui=confs[,3],xlab=expression(paste("Initial vale of ", mu)),ylab="Coverage probability")
  abline(h = 0.95, col = "red")
  
}


imu <- c(2,3)
itheta <- c(1,4)
isx <- c(3,1)
isy <- c(4,2)
int_theta <-   seq(0,18,2)

load("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git/data1001/Init_theta_vari_100.RData")
poOLD <- parms_org
sdOLD <- sd_org

load("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git/data1001/TMB_init_theta_varieret_100.RData")
poTMB <- parms_org
sdTMB <- sd_org


###### mean and sd of all parms ------------
#mu
par(mfrow=c(2,1),mai = c(0.8, 0.8, 0.1, 0.1))
plot(int_theta,apply((poOLD[,3,]),2,mean),pch=15,col="red",ylab = expression(paste('Mean of estimated vaue of ',mu)),xlab = expression(paste('Initial vale of ', theta)))
points(int_theta,apply((poTMB[,2,]),2,mean),col="blue",pch=8)
legend("top",c('Org','TMB'),pch=c(15,8),col=c('red','blue'))
abline(h =1)

plot(int_theta,apply((poOLD[,3,]),2,sd),pch=15, col="red",ylab = expression(paste('Sd of estimated vaue of ',mu)),xlab = expression(paste('Initial vale of ', theta)))
points(int_theta,apply((poTMB[,2,]),2,sd),col="blue",pch=8)
legend("top",c('Org','TMB'),pch=c(15,8),col=c('red','blue'))

#theta
plot(int_theta,apply((poOLD[,4,]),2,mean),pch=15,col="red",ylab = expression(paste('Mean of estimated vaue of ',theta)),xlab = expression(paste('Initial vale of ', theta)))
legend("bottomleft",c('Org','TMB'),pch=c(15,8),col=c('red','blue'))
points(int_theta,apply(exp(poTMB[,1,]),2,mean),col="blue",pch=8)
abline(h =10)

plot(int_theta,apply((poOLD[,4,]),2,sd),pch=15, col="red",ylab = expression(paste('Sd of estimated vaue of ',theta)),xlab = expression(paste('Initial vale of ', theta)))
legend("topleft",c('Org','TMB'),pch=c(15,8),col=c('red','blue'))
points(int_theta,apply(exp(poTMB[,1,]),2,sd),col="blue",pch=8)

#sigx
plot(int_theta,apply(exp(poOLD[,1,]),2,mean),pch=15,col="red",ylim=c(0.96,1.03),ylab = expression(paste('Mean of estimated vaue of ',sigma[x])),xlab = expression(paste('Initial vale of ', theta)))
legend("bottomleft",c('Org','TMB'),pch=c(15,8),col=c('red','blue'))
points(int_theta,apply(exp(poTMB[,3,]),2,mean),col="blue",pch=8)
abline(h =1)

plot(int_theta,apply(exp(poOLD[,1,]),2,sd),pch=15, col="red",ylab = expression(paste('Sd of estimated vaue of ',sigma[x])),xlab = expression(paste('Initial vale of ', theta)))
legend("topleft",c('Org','TMB'),pch=c(15,8),col=c('red','blue'))
points(int_theta,apply(exp(poTMB[,3,]),2,sd),col="blue",pch=8)

#sigy

plot(int_theta,apply(exp(poOLD[,2,]),2,mean),pch=15,col="red",ylim=c(0.009,0.0103),ylab = expression(paste('Mean of estimated vaue of ',bold(sigma[y]))),xlab = expression(paste('Initial vale of ', theta)))
legend("bottomleft",c('Org','TMB'),pch=c(15,8),col=c('red','blue'))
points(int_theta,apply(exp(poTMB[,4,]),2,mean),col="blue",pch=8)
abline(h =0.01)

plot(int_theta,apply(exp(poOLD[,2,]),2,sd),pch=15, col="red",ylim = c(0.0017,0.0028),ylab = expression(paste('Sd of estimated vaue of ',bold(sigma[y]))),xlab = expression(paste('Initial vale of ', theta)))
legend("topleft",c('Org','TMB'),pch=c(15,8),col=c('red','blue'))
points(int_theta,apply(exp(poTMB[,4,]),2,sd),col="blue",pch=8)



#### sigma fordeling 
#det er den samme for alle parametre jo.
sigma_stat <- matrix(data = NA, ncol = 3, nrow = 10)
s <- sdold
parameter <- 3
sigma_stat[,1] <- apply(s[,parameter,],2,mean)
sigma_stat[,2] <- apply(s[,parameter,],2,median)
sigma_stat[,3] <- apply(s[,parameter,],2,sd)
View(sigma_stat)


####### cov prob plots --------------


#looking at mu
par(mfrow=c(2,2),mai = c(0.8, 0.8, 0.5, 0.1))
CICI(poOLD,imu[2],sdOLD,int_theta,1,5)
title("Ctsmr",cex.main=0.8)
CICI(poTMB,imu[1],sdTMB,int_theta,1,5) 
xlab(expression(paste("Initial value of "), theta))
title("CtsmrTMB",cex.main=0.8)

mtext(expression(paste(mu, " when initial guess of ", theta, " is changed")),side=3,outer=TRUE, line=-1.5,cex.main=1.5)

#looking at theta
#par(mfrow=c(2,1))
CICI(poOLD,itheta[2],sdOLD,int_theta,10,5)
title("Ctsmr",cex.main=0.8)
CICI(poTMB,itheta[1],sdTMB,int_theta,log(10),5)
title("CtsmrTMB",cex.main=0.8)
mtext(expression(paste(theta, " when initial guess of ", theta, " is changed")),side=1,outer=TRUE, line=-15,cex.main=1.5)

#lookign at sig x
par(mfrow=c(2,2))
CICI(poOLD,isx[2],sdOLD,int_theta,log(1),2)
title("Ctsmr",cex.main=0.8)
CICI(poTMB,isx[1],sdTMB,int_theta,log(1),2)
title("CtsmrTMB",cex.main=0.8)
mtext(expression(paste(sigma[x], " when initial guess of ", theta, " is changed")),side=3,outer=TRUE, line=-1.5,cex.main=1.5)

#looking at sigy
#par(mfrow=c(1,2))
CICI(poOLD,isy[2],sdOLD,int_theta,log(1e-2),1)
title("Ctsmr",cex.main=0.8)
CICI(poTMB,isy[1],sdTMB,int_theta,log(1e-2),1)
title("CtsmrTMB",cex.main=0.8)
mtext(expression(paste(sigma[y], " when initial guess of ", theta, " is changed")),side=1,outer=TRUE, line=-15,cex.main=1.5)

### count nan/na ----

