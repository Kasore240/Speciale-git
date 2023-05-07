OU_ctsmrTMB <- function(init_pars,init_lb,init_ub, sx0){
  
  obj = ctsmrTMB$new()
  
  obj$set_modelname("ornstein_uhlenbeck")
  
  # Add system equations
  obj$add_systems(
    dx ~ theta * (mu-x) * dt + sigma_x*dw
  )
  
  # Add observation equations
  obj$add_observations(
    y ~ x
  )
  
  # Set observation equation variances
  obj$add_observation_variances(y ~ sigma_y^2)
  
  # Specify algebraic relations   -- er de her mon overflødige?
  obj$add_algebraics(
    theta ~ exp(logtheta),
    sigma_x ~ exp(logsigma_x),
    sigma_y ~ exp(logsigma_y)
  )
  
  # Specify parameter initial values and lower/upper bounds in estimation
  obj$add_parameters(
    logtheta ~ log(c(init = init_pars[1], lower=init_lb[1], upper=init_ub[1])),
    mu ~ c(init=init_pars[2], lower=init_lb[1], upper=init_ub[2]),
    logsigma_x ~ log(c(init= init_pars[3], lower=init_lb[3], upper=init_ub[3])),
    logsigma_y ~ log(c(init=init_pars[4], lower=init_lb[4], upper=init_ub[4]))
  )
  # Set initial state mean and covariance
  obj$set_initial_state(init_pars[5], sx0*diag(1))
  
  return(obj)}