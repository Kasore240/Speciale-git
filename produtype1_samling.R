library(ctsmrTMB)
library(ggplot2)
##### --------------------------------------------------######
# 
obj = ctsmrTMB$new()

# Modelname
obj$set_modelname("ornstein_uhlenbeck")

# Set location path for C++ file generation
# obj$set_path("~/")

# Add system equations
obj$add_systems(
  dx ~ theta * (mu-x) * dt + sigma_x*dw
)

# Add observation equations
obj$add_observations(
  y ~ x
)

# Set observation equation variances
obj$add_observation_variances(
  y ~ sigma_y^2
)

# Specify algebraic relations
obj$add_algebraics(
  theta ~ exp(logtheta),
  sigma_x ~ exp(logsigma_x),
  sigma_y ~ exp(logsigma_y)
)

# Perform a lamperti transformation (for state dependent diffusion)
# obj$set_lamperti("log")

# Specify inputs
# obj$add_inputs(input1, input2)

# Specify parameter initial values and lower/upper bounds in estimation
obj$add_parameters(
  logtheta ~ log(c(init = 1, lower=1e-5, upper=50)),
  mu ~ c(init=1.5, lower=0, upper=5),
  logsigma_x ~ log(c(init= 1e-1, lower=1e-10, upper=10)),
  logsigma_y ~ log(c(init=1e-1, lower=1e-10, upper=10))
)
N.sim <- 1000

# Simulate data using Euler Maruyama
set.seed(10)
pars = c(theta=10, mu=1, sigma_x=1, sigma_y=1e-2)
# 
dt.sim = 1e-3
t.sim = seq(0,1,by=dt.sim)
dw = rnorm(length(t.sim)-1,sd=sqrt(dt.sim))
x = 3
for(i in 1:N.sim) {
  x[i+1] = x[i] + pars[1]*(pars[2]-x[i])*dt.sim + pars[3]*dw[i]
}

# Extract observations and add noise
dt.obs = 1e-2
t.obs = seq(0,1,by=dt.obs)
y = x[t.sim %in% t.obs] + sqrt(pars[4]) * rnorm(length(t.obs))

# Set initial state mean and covariance
obj$set_initial_state(x[1], 1e-1*diag(1))

# Create data
.data = data.frame(
  t = t.obs,
  y = y
)

# Carry out estimation using extended kalman filter method
fit <- obj$estimate(.data, method="ekf", use.hessian=T)

# Check parameter estimates against truth
pars2real = function(x) c(exp(x[1]),x[2],exp(x[3]),exp(x[4]))
sum_phillip <- cbind( pars2real(fit$par.fixed), pars )
sum_phillip
# plot one-step predictions, simulated states and observations
t.est = fit$states$mean$prior$t
x.mean = fit$states$mean$prior$x
x.sd = fit$states$sd$prior$x
ggplot() +
  geom_ribbon(aes(x=t.est, ymin=x.mean-2*x.sd, ymax=x.mean+2*x.sd),fill="grey", alpha=0.9) +
  geom_line(aes(x=t.sim,y=x)) + 
  geom_line(aes(x=t.est, x.mean),col="blue") +
  geom_point(aes(x=t.obs,y=y),col="red",size=2) +
  theme_minimal()


# Check one-step-ahead residuals
plot(fit, use.ggplot=F)

# Use prediction function to get 10-step-ahead state predictions
pred = predict(fit, n.step.ahead=10)

##### --------------------------------------------------######
#----------------------------------------------------------------
# Init by deleting all variables and functions
#rm(list=ls())
# Set the working directory
#setwd("C:/Users/bruger/OneDrive/Skrivebord/Speciale")

# Use the ctsmr package
library(ctsmr)
#----------------------------------------------------------------

#----------------------------------------------------------------
# Generate a new object of class ctsm.
model <- ctsm$new()

# Now we can add a system equation (and thereby also a state variable) with the addSystem() function
#model$addSystem(dTi ~ ( 1/(Ci*Ria)*(Ta-Ti) + gA/Ci*Ps + 1/Ci*Ph )*dt + exp(p11)*dw1)

model$addSystem(
  dx ~ theta * (mu-x) * dt + exp(log_sigma_x)*dw
)


# Set the observation equation: Ti is the state, yTi is the measured output
#model$addObs(yTi ~ Ti)
model$addObs(y ~ x)

# Set the variance of the measurement error
#model$setVariance(yTi ~ exp(e11))
model$setVariance(y ~ exp(logsigma2_y))

# Set the initial value of the value of the state at the start time (values where
# the estimation (i.e. optimization of the likelihood) starts) and also the lower
# and upper bound, which must contain the parameter value
#model$setParameter(  Ti = c(init=15  ,lb=0     ,ub=25 ) )
model$setParameter(
  theta = c(init = 1, lb=1e-5, ub=50),
  mu = c(init=1.5, lb=0, ub=5),
  log_sigma_x = log(c(init= 1e-1, lb=1e-10, ub=10)),
  logsigma2_y = log(c(init=1e-1, lb=1e-10, ub=10))
)

model$setParameter(x0 = c(init=10,lb=1,ub=100))

# Run the parameter estimation
fit <- model$estimate(.data)

# See the summary of the estimation
sum_ctsr <- summary(fit)
sum_ctsr$coefficients[2:3,'Estimate'] <- exp(sum_ctsr$coefficients[2:3,'Estimate'])
sum_ctsr$coefficients[3,'Estimate'] <- sqrt(sum_ctsr$coefficients[3,'Estimate'])
#----------------------------------------------------------------
sum_ctsr$coefficients[,'Estimate'] 
sum_phillip
#----------------------------------------------------------------
# If any of the parameter estimates are close to the lower or upper bound then "dPen/dPar" is
# significant compered to "dF/dPar"
summary(fit, extended = TRUE)
#----------------------------------------------------------------


#----------------------------------------------------------------
# Calculate the one-step predictions of the state (i.e. the residuals)
tmp <- predict(fit)[[1]]
str(tmp)

# Calculate the residuals and put them with the data in a data.frame X
.data$residuals <- .data$y - tmp$output$pred$y

# Plot the auto-correlation function and cumulated periodogram in a new window
par(mfrow=c(1,3))
# The blue lines indicates the 95% confidence interval, meaning that if it is
#  white noise, then approximately 19 out of 20 lag correlations will be inside.
acf(.data$residuals, lag.max=8*24)
# The periodogram is the estimated energy spectrum in the signal
spec.pgram(.data$residuals)
# The cumulated periodogram 
cpgram(.data$residuals)
#----------------------------------------------------------------
plot(.data$residuals,type='l')

