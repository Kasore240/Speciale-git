#load("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git/fordeling startgæt/fordeling_bigstep7000_cluster.Rdata")
load("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git/fordeling startgæt/fordeling_smallstep7000_cluster.Rdata")
load("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git/fordeling startgæt/fordeling_bigstep7000_cluster.Rdata")
load("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git/fordeling startgæt/fordeling_smallstep_teo7000_cluster.Rdata")
load("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git/fordeling startgæt/fordeling_bigstep_teo7000_cluster.Rdata")
#sig X
par(mfrow=c(2,1))
hist(exp(po[,1]),breaks = 40,main='ctsmr',cex.main=0.8, xlab=expression(paste("Estimated ", widetilde(sigma)[x])))
hist(exp(pt[,3]),breaks=40,main = "ctsmrTMB",cex.main=0.8,xlab=expression(paste("Estimated ", widetilde(sigma)[x])))
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


i <- 2
c(sd((pt[,i])),mean((pt[,i])),sd((spt[,i])),mean((spt[,i])))

########## coverage ###############
p <- 1e-2
i <- 2
po <- po[so[,2]<2,]
so <- so[so[,2]<2,]
low =exp(c(po[,i]) - (2*c(so[,i])))
up =exp(c(po[,i]) + (2*c(so[,i])))
x <- sum(p< up & p> low)
b <- binom.test(x,length(po[,1]),0.95,conf.level = 0.95)
c(b$conf.int[1],b$estimate,b$conf.int[2])

load("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git/fordeling startgæt/fordeling_smallstep_teo7000_cluster.Rdata")
load("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git/fordeling startgæt/fordeling_bigstep_teo7000_cluster.Rdata")
colMeans((st))
apply((st),2,median)
apply((st),2,sd)

c( sum(so[,1]>2) , sum(st[,4]>2) )
sum(so[,2]>2 & st[,4]>2)
