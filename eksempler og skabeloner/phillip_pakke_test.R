library(ctsmrTMB)
library(ggplot2)

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
cbind( pars2real(fit$par.fixed), pars )

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
