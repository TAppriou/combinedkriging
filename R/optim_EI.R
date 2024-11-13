optim_EI <- function(model, upper, lower, full.exploitation=FALSE, combination=TRUE, optimizer="cmaes") {

  if (optimizer!="cmaes") {
    stop("Only cmaes implemented for this version")
  }

  # Parameters of the CMA-ES
  d <- model$d
  lambda <- min(32*d, 32*50) #Population size
  maxit <- 50
  parinit <- runif(d, min=lower, max=upper)

  # Optimization
  opt <- cma_es_bis(par=parinit, fn=compute_EI, model=model, full.exploitation=full.exploitation, combination=combination,
                    upper=upper, lower=lower,
                    control = list(maxit=maxit, lambda=lambda, vectorized=TRUE))


  return(list(par=opt$par, value=opt$value))

}
