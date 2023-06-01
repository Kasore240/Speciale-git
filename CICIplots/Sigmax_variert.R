library("plotrix")
CICI <- function(po,index,sd,xscale,p,v){
  confs <- matrix(data = NA, nrow = 10, ncol = 3)
  for (k in 1:10) {
    if (length(p) > 1){ tp <- p[k]}
    else{ tp <- p}
    
    pol <- po[sd[,index,k]<v,index,k]
    sdl <- sd[sd[,index,k]<v,index,k]
    low =(c(pol) - (qt(0.975,length( pol ))*c(sdl)))
    up =(c(pol) + (qt(0.975,length( pol ))*c(sdl)))
    x <- sum(tp< up & tp> low,na.rm=T)
    
    b <- binom.test(x,length(low),0.95,conf.level = 0.95)
    confs[k,1] <- x/length(low)
    confs[k,2:3] <- b$conf.int
  }
  
  plotCI(x=xscale, y= confs[,1],li=confs[,2],ui=confs[,3],xlab=expression(paste('Size of ', sigma[x])),ylab="Coverage probability")
  abline(h = 0.95, col = "red")
  
}

imu <- c(2,3)
itheta <- c(1,4)
isx <- c(3,1)
isy <- c(4,2)
int_sig_x <- seq(1e-1,5,0.49)
int_sig_x <- int_sig_x[1:10]
load("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git/small500/sigmax_varieret_100.RData")
poOLD1 <- parms_org
sdOLD1 <- sd_org

poTMB1 <- parms_tmb
sdTMB1 <- sd_tmb

apply(is.na(sdOLD[,1,]),2,sum)
apply(is.na(sdTMB[,1,]),2,sum)

apply((sdOLD1[,2,]>2),2,sum,na.rm=TRUE)
apply((sdTMB1[,4,]>2),2,sum,na.rm=TRUE)

for (i in 1:10){
  for (j in 1:4){
  sdOLD[sdOLD[,2,i]>2,j,i] <- NA
  sdTMB[sdTMB[,4,i]>2,j,i] <- NA
  poOLD[sdOLD[,2,i]>2,j,i] <- NA
  poTMB[sdTMB[,4,i]>2,j,i] <- NA
}}

###### mean and sd of all parms ------------
#mu
par(mfrow=c(2,1),mai = c(0.8, 0.8, 0.1, 0.1))
plot(int_sig_x ,apply((poOLD[,3,]),2,mean,na.rm=TRUE),pch=15,col="red",ylim = c(0.95,1.1),ylab = expression(paste('Mean of estimated vaue of ',mu)),xlab = expression(paste('Size of ', sigma[x])))
points(int_sig_x ,apply((poTMB[,2,]),2,mean,na.rm=TRUE),col="blue",pch=8)
legend("topleft",c('Ctsmr','CtsmrTMB'),pch=c(15,8),col=c('red','blue'))
abline(h =1)

plot(int_sig_x ,apply((poOLD[,3,]),2,sd,na.rm=TRUE),pch=15, col="red",ylab = expression(paste('Sd of estimated vaue of ',mu)),xlab = expression(paste('Size of ', sigma[x])))
points(int_sig_x ,apply((poTMB[,2,]),2,sd,na.rm=TRUE),col="blue",pch=8)
legend("topleft",c('Ctsmr','CtsmrTMB'),pch=c(15,8),col=c('red','blue'))

#theta
plot(int_sig_x ,apply((poOLD[,4,]),2,mean,na.rm=TRUE),pch=15,col="red",ylim = c(10,14), ylab = expression(paste('Mean of estimated vaue of ',theta)),xlab = expression(paste('Size of ', sigma[x])))
points(int_sig_x ,apply(exp(poTMB[,1,]),2,mean,na.rm=TRUE),col="blue",pch=8)
legend("topleft",c('Ctsmr','CtsmrTMB'),pch=c(15,8),col=c('red','blue'))
abline(h =10)

plot(int_sig_x ,apply((poOLD[,4,]),2,sd,na.rm=TRUE),pch=15, col="red",ylim=c(0,6.5),ylab = expression(paste('Sd of estimated vaue of ',theta)),xlab = expression(paste('Size of ', sigma[x])))
legend("topleft",c('Ctsmr','CtsmrTMB'),pch=c(15,8),col=c('red','blue'))
points(int_sig_x ,apply(exp(poTMB[,1,]),2,sd,na.rm=TRUE),col="blue",pch=8)

#sigx
plot(int_sig_x ,apply(exp(poOLD[,1,]),2,mean,na.rm=TRUE),pch=15,col="red",ylim=c(0,5),ylab = expression(paste('Mean of estimated vaue of ',sigma[x])),xlab = expression(paste('Size of ', sigma[x])))
legend("topleft",c('Ctsmr','CtsmrTMB'),pch=c(15,8),col=c('red','blue'))
points(int_sig_x ,apply(exp(poTMB[,3,]),2,mean,na.rm=TRUE),col="blue",pch=8)
abline(0,1)

plot(int_sig_x ,apply(exp(poOLD[,1,]),2,sd,na.rm=TRUE),pch=15, col="red",ylim=c(0,0.2),ylab = expression(paste('Sd of estimated vaue of ',sigma[x])),xlab = expression(paste('Size of ', sigma[x])))
legend("left",c('Ctsmr','CtsmrTMB'),pch=c(15,8),col=c('red','blue'))
points(int_sig_x ,apply(exp(poTMB[,3,]),2,sd,na.rm=TRUE),col="blue",pch=8)

#sigy

plot(int_sig_x ,apply(exp(poOLD[,2,]),2,mean,na.rm=TRUE),pch=15,col="red",ylim=c(0,0.02),ylab = expression(paste('Mean of estimated vaue of ',bold(sigma[y]))),xlab = expression(paste('Size of ', sigma[x])))
legend("topleft",c('Ctsmr','CtsmrTMB'),pch=c(15,8),col=c('red','blue'))
points(int_sig_x ,apply(exp(poTMB[,4,]),2,mean,na.rm=TRUE),col="blue",pch=8)
abline(h =0.01)

plot(int_sig_x ,apply(exp(poOLD[,2,]),2,sd,na.rm=TRUE),pch=15, col="red",ylim=c(0,0.05),ylab = expression(paste('Sd of estimated vaue of ',bold(sigma[y]))),xlab = expression(paste('Size of ', sigma[x])))
legend("topleft",c('Ctsmr','CtsmrTMB'),pch=c(15,8),col=c('red','blue'))
points(int_sig_x ,apply(exp(poTMB[,4,]),2,sd,na.rm=TRUE),col="blue",pch=8)



#### sigma fordeling  ##########
#det er den samme for alle parametre jo.
sigma_stat <- matrix(data = NA, ncol = 4, nrow = 10)
s <- sdOLD   #old 2, tmb 4
parameter <- 1
sigma_stat[,1] <- apply(s[,parameter,],2,mean,na.rm=TRUE)
sigma_stat[,2] <- apply(s[,parameter,],2,median,na.rm=TRUE)
sigma_stat[,3] <- apply(s[,parameter,],2,sd,na.rm=TRUE)
sigma_stat[,4] <- apply(s[,parameter,],2,mean,na.rm=TRUE) -apply(s[,parameter,],2,median,na.rm=TRUE)
View(sigma_stat)


######### cov plot ################



#looking at mu
par(mfrow=c(2,2),mai = c(0.8, 0.8, 0.5, 0.1))
CICI(poOLD,imu[2],sdOLD,int_sig_x,1,5)
title("Ctsmr",cex.main=0.8)
CICI(poTMB,imu[1],sdTMB,int_sig_x,1,5) 
title("CtsmrTMB",cex.main=0.8)

mtext(expression(paste(mu, " when size of ", sigma[x], " is changed")),side=3,outer=TRUE, line=-1.5,cex.main=1.5)

#looking at theta
#par(mfrow=c(2,1))
CICI(poOLD,itheta[2],sdOLD,int_sig_x,10,5)
title("Ctsmr",cex.main=0.8)
CICI(poTMB,itheta[1],sdTMB,int_sig_x,log(10),5)
title("CtsmrTMB",cex.main=0.8)
mtext(expression(paste(theta, " when size of ", sigma[x], " is changed")),side=1,outer=TRUE, line=-18,cex.main=1.5)

#lookign at sig x
par(mfrow=c(2,2))
CICI(poOLD,isx[2],sdOLD,int_sig_x,log(int_sig_x),2)
title("Ctsmr",cex.main=0.8)
CICI(poTMB,isx[1],sdTMB,int_sig_x,log(int_sig_x),2)
title("CtsmrTMB",cex.main=0.8)
mtext(expression(paste(sigma[x], " when size of ", sigma[x], " is changed")),side=3,outer=TRUE, line=-1.5,cex.main=1.5)

#looking at sigy
#par(mfrow=c(1,2))
CICI(poOLD,isy[2],sdOLD,int_sig_x,log(1e-2),1)
title("Ctsmr",cex.main=0.8)
CICI(poTMB,isy[1],sdTMB,int_sig_x,log(1e-2),1)
title("CtsmrTMB",cex.main=0.8)
mtext(expression(paste(sigma[y], " when size of ", sigma[x], " is changed")),side=1,outer=TRUE, line=-18,cex.main=1.5)




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
