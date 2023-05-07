sim_OU_EM <- function(N.sim, dt.sim, dt.obs, pars,x0){
  
  t.sim = seq(0,1,by=dt.sim)
  dw = rnorm(length(t.sim)-1,sd=sqrt(dt.sim))
  x = x0
  for(i in 1:(length(t.sim)-1)) {
    x[i+1] = x[i] + pars[1]*(pars[2]-x[i])*dt.sim + pars[3]*dw[i]*sqrt(x[i])
  }
  
  # Extract observations and add noise
  t.obs = seq(0,1,by=dt.obs)
  y = x[seq(1,length(x),dt.obs/dt.sim)]+ pars[4] * rnorm(length(t.obs))
  
  
  
  # Create data
  .data = data.frame(
    t = t.sim[seq(1,length(x),dt.obs/dt.sim)],
    y = y
  )
  l1 <- list(.data=.data,x=x)
  return(l1)}