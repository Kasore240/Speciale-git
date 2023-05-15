#load("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git/fordeling.Rdata")
load("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git/fordeling_teoretisk.Rdata")
#sig X
par(mfrow=c(2,1))
hist(exp(po[,1]),breaks = 40,main='ctsmr',cex.main=0.8, xlab=expression(paste("Estimated ", widetilde(sigma)[x])))
hist(exp(pt[,3]),breaks=40, main = "ctsmrTMB",cex.main=0.8,xlab=expression(paste("Estimated ", widetilde(sigma)[x])))
mtext(expression(paste('Distribution of ', widetilde(sigma)[x])),side=3,outer=TRUE, line=-1.5,cex.main=1.4)

#sig y
par(mfrow=c(2,1))
hist(exp(po[,2]),breaks=50,main='ctsmr',cex.main=0.8, xlab=expression(paste("Estimated ", widetilde(sigma)[y])))
hist(exp(pt[,4]),breaks=50,main='ctsmrTMB',cex.main=0.8, xlab=expression(paste("Estimated ", widetilde(sigma)[y])))
mtext(expression(paste('Distribution of ', widetilde(sigma)[y])),side=3,outer=TRUE, line=-1.5,cex.main=1.4)

#mu
par(mfrow=c(2,1))
hist((po[,3]),breaks=40,main='ctsmr',cex.main=0.8, xlab=expression(paste("Estimated ", widetilde(mu))))
hist((pt[,2]),breaks=40,main='ctsmrTMB',cex.main=0.8, xlab=expression(paste("Estimated ", widetilde(mu))))
mtext(expression(paste('Distribution of ', widetilde(mu))),side=3,outer=TRUE, line=-1.5,cex.main=1.4)
#theta
par(mfrow=c(2,1))
hist((po[,4]),breaks=40,main='ctsmr',cex.main=0.8, xlab=expression(paste("Estimated ", widetilde(theta))))
hist(exp(pt[,1]),breaks=40,main='ctsmrTMB',cex.main=0.8, xlab=expression(paste("Estimated ", widetilde(theta))))
mtext(expression(paste('Distribution of ', widetilde(theta))),side=3,outer=TRUE, line=-1.5,cex.main=1.4)
