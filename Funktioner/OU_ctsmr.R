OU_ctsmr <- function(init_pars,init_lb,init_ub){
  
  model <- ctsm$new()
  
  # Add system equations
  model$addSystem( dx ~ theta * (mu-x) * dt + exp(log_sigma_x)*dw)
  
  # Add observation equations
  model$addObs(y ~ x)
  
  # Set observation equation variances
  model$setVariance(y ~ (exp(logsigma2_y))^2)
  
  model$setParameter(x0 = c(init=init_pars[5],lb=init_lb[5],ub=init_ub[5]))
  
  # Specify parameter initial values and lower/upper bounds in estimation
  model$setParameter(
    logsigma2_y = log(c(init=init_pars[4], lb=init_lb[4], ub=init_ub[4])),
    log_sigma_x = log(c(init= init_pars[3], lb=init_lb[3], ub=init_ub[3])),
    mu = c(init=init_pars[2], lb=init_lb[2], ub=init_ub[2]),
    theta = c(init = init_pars[1], lb=init_lb[1], ub=init_ub[1])
  )
  return(model)
}