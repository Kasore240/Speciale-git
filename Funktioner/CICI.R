CICI <- function(po,index,sd,xscale,p){
  confs <- matrix(data = NA, nrow = 10, ncol = 3)
  for (k in 1:10) {
    low =(c(po[,index,k]) - (2*c(sd[,index,k])))
    up =(c(po[,index,k]) + (2*c(sd[,index,k])))
    x <- sum(p[k]< up & p[k]> low,na.rm=T)
    a <- 100 - sum(is.na(sd[,index,k]))
    b <- binom.test(x,a,0.95,conf.level = 0.95)
    
    confs[k,1] <- x/a
    confs[k,2:3] <- b$conf.int
  }
  
  plotCI(x=xscale, y= confs[,1],li=confs[,2],ui=confs[,3])
  abline(h = 0.95, col = "red")
  
}