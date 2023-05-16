# setwd("~/Documents/imm/phd/2015_summer_school_DTU/rj/presentation/rj/likelihood")

# Load a few required packages
#install.packages("numDeriv")
#install.packages("animation")

library(numDeriv)
library(animation)



# Specify the population parameters
mu <- 150
sigma <- 15
# Number of samples
n <- 100
# Randomly sample n observations the normal distribution
y <- rnorm(n = n, mean = mu, sd = sigma)

# See the sample
hist(y, prob=TRUE)
xseq <- seq(50,250,len=1000)
lines(xseq, dnorm(xseq, mean = mu, sd = sigma))

# The negative log-likelihood as a function of the parameters theta
neg_loglike <- function(theta) {
    mu <- theta[1]
    var <- exp(theta[2])
    #
    ll <- n*log(sqrt(var)) + 1/2 * sum((y - mu)^2/var)
    # Another way:   ll <- -sum(dnorm(y, mean=mu, sd=sqrt(var),log=TRUE))
    return(ll)
}


# Run an optimization of the negative log-likelihood function

# Optimise the function starting at (0,0).

# nlminb outputs the current point at every iteration. The output is captured and

# saved in a temporary file.
ss <- capture.output(opt <- nlminb(start=c(0,0), objective=neg_loglike, control=list(trace=1)),
                     file = out.file <- tempfile())
                                        # Okay.. Load it here
trace <- read.table(out.file, stringsAsFactors=FALSE)[,-1]
colnames(trace) <- c("ll","mu","lvar")
trace$ll <- as.double(substr(trace$ll, 1L, 9L))

ymax <- dnorm(trace$mu, mean=trace$mu, sd = sqrt(exp(trace$lvar)))[nrow(trace)]

# Animation running the optimisation
oopt <- ani.options(interval = .1)

px <- seq(-200,400,length.out=400)

opar <- par(mfrow=c(1,2), oma=c(3,3,4,0), mar=c(4,3,1,1), las=1, cex=0.8)

for (i in seq_len(nrow(trace))) {
    pdf <- dnorm(x=px, mean=trace$mu[i], sd=sqrt(exp(trace$lvar[i])))
    #
    plot(px, pdf, type="l", xlab = "Height", ylab = "PDF", ylim = c(0, ymax))
    lines(px, dnorm(x=px, mean=mu, sd=sigma), col="blue")
    rug(y, col="red", lwd=2)
    #  
    plot(seq_len(i), trace$ll[1:i], 
         xlim = c(1,nrow(trace)), ylim = range(trace$ll), type="p", 
         xlab = "Iteration", log="y", pch=16)
    #  
    #  
    title(main = bquote(atop(list(Iteration == .(i), 
                                  -log(L) == .(formatC(trace$ll[i], digits=1, format="f"))
                                  ), 
                             list(mu    == .(formatC( trace$mu[i], digits=1, format="f")),
                                  sigma == .(formatC( sqrt(exp(trace$lvar[i])), digits=1, format="f"))
                                  )
                             )
                        ), outer = TRUE)
    #
    ani.pause(0.4)
}

# Restore the plotting pars.
par(opar)


# Contour plot of the neg_loglikelihood function of the parameter space.

mu_seq <- seq(-20,300, length.out=100)
lsd <- seq(-3, 12, length.out=100)
#
grid <- as.matrix(expand.grid(mu = mu_seq, lsd = lsd))
#
llc <- numeric(nrow(grid))
#
for (i in seq_len(nrow(grid))) {
    llc[i] <- neg_loglike(grid[i,])
}
#
for (i in seq_len(nrow(trace))) {
    contour(x = mu_seq,
            y = lsd,
            z = matrix(log(llc),
                       ncol=length(mu_seq)), 
            levels = c(seq(log(min(llc)), 8, length.out=20),
                       seq(8,20,length.out=30)),
            xlab = bquote(mu_seq),
            ylab = bquote(log(sigma^2)))
    #  
    points(mu,log(sigma^2), col="blue", pch=16)
    #  
    lines(trace$mu[seq_len(i)], trace$lvar[seq_len(i)], col="red",lwd=2)
    points(trace$mu[seq_len(i)], trace$lvar[seq_len(i)], col="red",pch=16)
    #  
    title(main = bquote(atop(list(Iteration == .(i), 
                                  -log(L) == .(formatC(trace$ll[i], digits=1, format="f"))
                                  ), 
                             list(mu    == .(formatC( trace$mu[i], digits=1, format="f")),
                                  sigma == .(formatC( sqrt(exp(trace$lvar[i])), digits=1, format="f"))
                                  )
                             )
                        ))
    #  
    ani.pause(0.2)
}



#----------------------------------------------------------------
# The estimated values

# The population mean
mu
# The estimated mean
opt$par[1]

# The population standard deviation
sigma
# The estimated mean
sqrt(exp(opt$par[2]))


#----------------------------------------------------------------
# Use linear regression and compare
fit <- lm(y ~ 1)
summary(fit)
sd(fit$residuals)


#----------------------------------------------------------------
# What about the estimated uncertainty of the parameters?

# I.e. the standard error ("Std. Error" in)
summary(fit)

# Approximate the Hessian
H <- hessian(neg_loglike, opt$par)

# Fisher's information matrix
I <- qr.solve(H)

# Covariance matrix of the estimates
cov2cor(I)

# Standard error (Standard deviation estimate of the parameters)
sqrt(diag(I))
