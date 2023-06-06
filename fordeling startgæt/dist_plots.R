#load("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git/fordeling startgæt/fordeling_bigstep7000_cluster.Rdata")
load("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git/true data/fordeling_smallstep_teo7000_cluster.Rdata")
load("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git/true data/fordeling_smallstep7000_cluster.Rdata")
load("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git/true data/fordeling_bigstep_teo7000_cluster.Rdata")
load("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git/true data/fordeling_bigstep7000_cluster.Rdata")
#sig X

b <- min(c(exp(po[,1]),exp(pt[,3]))) - 0.001 # Set the minimum for the breakpoints
e <- max(c(exp(po[,1]),exp(pt[,3])))+0.5 # Set the maximum for the breakpoints
ax <- seq(e, n = 15)
c1 <- rgb(173,216,230,max = 255, alpha = 125, names = "blue")
c2 <- rgb(255,192,203, max = 255, alpha = 125, names = "red")

par(mfrow=c(2,1))
ho <-hist(exp(po[,1]),breaks=45,plot=FALSE)
ht <- hist(exp(pt[,3]),breaks=45,plot=FALSE)
plot(ho, col = c1,xlab=expression(paste("Estimated ", sigma[x])),cex.lab=1.2,main = expression(paste('Distribution of estimated ', sigma[x])),ylim = c(0,1000)) 
plot(ht, col = c2, add = TRUE)
legend("topleft",legend = c("Ctsmr","CtsmrTMB"),pch=c(15,15),col = c("blue","red"))


#sig y

ho <-hist(exp(po[,2]),breaks=45,plot=FALSE)
ht <- hist(exp(pt[,4]),breaks=45,plot=FALSE)
plot(ht, col = c2,xlab=expression(paste("Estimated ", sigma[y])),cex.lab=1.2,main = expression(paste('Distribution of estimated ', sigma[y])),ylim = c(0,1100)) 
plot(ho, col = c1, add = TRUE)
legend("topleft",legend = c("Ctsmr","CtsmrTMB"),pch=c(15,15),col = c("blue","red"))

#mu

ho <-hist((po[,3]),breaks=40,plot=FALSE)
ht <- hist((pt[,2]),breaks=40,plot=FALSE)
plot(ht, col = c2,xlab=expression(paste("Estimated ", mu)),cex.lab=1.2,main = expression(paste('Distribution of estimated ', mu)),ylim = c(0,800)) 
plot(ho, col = c1, add = TRUE)
legend("topleft",legend = c("Ctsmr","CtsmrTMB"),pch=c(15,15),col = c("blue","red"))
#theta

ho <-hist((po[,4]),breaks=40,plot=FALSE)
ht <- hist(exp(pt[,1]),breaks=40,plot=FALSE)
plot(ht, col = c2,xlab=expression(paste("Estimated ", theta)),cex.lab=1.2,main = expression(paste('Distribution of estimated ', theta)),ylim = c(0,900)) 
plot(ho, col = c1, add = TRUE)
legend("topleft",legend = c("Ctsmr","CtsmrTMB"),pch=c(15,15),col = c("blue","red"))


load("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git/true data/fordeling_smallstep_teo7000_cluster.Rdata")
load("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git/true data/fordeling_bigstep_teo7000_cluster.Rdata")

i <- 3
j <-2

c(sd((po[,i])),mean((po[,i])),sd((pt[,j])),mean((pt[,j])))

########## coverage ###############
p <- log(1e-2)
i <- 4

load("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git/true data/fordeling_bigstep7000_cluster.Rdata")
low =(c(pt[,i]) - (qt(0.975,length(st[,i]))*c(st[,i])))
up =(c(pt[,i]) + (qt(0.975,length(st[,i]))*c(st[,i])))
x <- sum(p< up & p> low,na.rm=TRUE)
b <- binom.test(x,length(st[,i]),0.95,conf.level = 0.95)
big <- c(b$conf.int[1],b$estimate,b$conf.int[2])

load("C:/Users/bruger/OneDrive/Skrivebord/Speciale/Speciale-git/true data/fordeling_smallstep7000_cluster.Rdata")
low =(c(pt[,i]) - (qt(0.975,length(st[,i]))*c(st[,i])))
up =(c(pt[,i]) + (qt(0.975,length(st[,i]))*c(st[,i])))
x <- sum(p< up & p> low,na.rm=TRUE)
b <- binom.test(x,length(st[,i]),0.95,conf.level = 0.95)
small <- c(b$conf.int[1],b$estimate,b$conf.int[2])

big
small





colMeans((st))
apply((st),2,median)
apply((st),2,sd)

c( sum(so[,1]>2) , sum(st[,3]>2,na.rm=TRUE) )
sum(so[,2]>2 & st[,4]>2)

mask <- st[,3]>2
mask[is.na(mask)] <- FALSE

med(st[!mask,3])
sigma_stat <- matrix(data = NA, ncol = 4, nrow = 4)
s <- st #old 2, tmb 4
parameter <- 1
sigma_stat[,1] <- apply(s,2,mean,na.rm=TRUE)
sigma_stat[,2] <- apply(s,2,median,na.rm=TRUE)
sigma_stat[,3] <- apply(s,2,sd,na.rm=TRUE)
sigma_stat[,4] <- apply(s,2,mean,na.rm=TRUE) -apply(s,2,median,na.rm=TRUE)
View(sigma_stat)

apply(so,2,sd,na.rm=TRUE)
