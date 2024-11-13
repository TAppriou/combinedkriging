krg_MLE <- function(X_train, Y_train, bornes_opt, kernel="matern5_2", theta_init=NULL, estim_mean=TRUE, nugget=1e-6) {

  if (estim_mean==TRUE) {
    trend <- NULL
  } else {
    trend <- 0
  }

  d <- dim(X_train)[2]

  # Kriging model
  krg_model <- DiceKriging::km(design=X_train, response=Y_train ,covtype=kernel, lower=rep(bornes_opt[1],d), upper=rep(bornes_opt[2],d),
                  estim.method="MLE", coef.trend=trend, parinit=theta_init, multistart=1, control=list(trace=FALSE, maxit=300), nugget=nugget)


  # Pre-compute useful quantities
  res <- list()

  res$X_train <- krg_model@X
  res$Y_train <- krg_model@y
  res$estim_mean <- estim_mean
  res$n_train <- dim(res$X_train)[1]
  res$d <- dim(res$X_train)[2]
  res$theta <- coef(krg_model,"range")
  res$kernel <- kernel

  # Covariance structure (for grad_EI)
  res$CovStruct <- krg_model@covariance

  # Cholesky of covariance matrix
  res$CholK <- krg_model@T

  # Mean of the GP
  res$mu <- krg_model@trend.coef

  # Cholesky times ones
  res$CholK_1 <- forwardsolve(t(res$CholK),matrix(1,nrow=res$n_train,ncol=1))

  # Cholesky times prediction
  res$CholK_Y <- krg_model@z

  # Variance of the GP
  res$sigma2 <- krg_model@covariance@sd2


  return(res)

}
