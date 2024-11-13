pred_krg <- function(X_test, model, fullCov=FALSE) {

  res <- list()

  if (is.matrix(X_test)) {
    n_test <- dim(X_test)[1]
  } else {
    n_test <- 1
  }

  # Covariances with prediction point
  K_X_Xtest <- model$sigma2 * nestedKriging::getCrossCorrMatrix(model$X_train, matrix(X_test, nrow=n_test, ncol=model$d), model$theta, model$kernel)
  CholK_K <- forwardsolve(t(model$CholK), K_X_Xtest)

  # Mean
  Y_pred <- matrix(model$mu, nrow=n_test, ncol=1) + crossprod(CholK_K, model$CholK_Y)
  res$mean <- as.vector(Y_pred)

  if (fullCov == FALSE) {
    # Variance
    var_pred <- rep(model$sigma2,length=n_test) - colSums(CholK_K^2)
    if (model$estim_mean==TRUE) {
      var_pred <- var_pred + (rep(1,length=n_test) - crossprod(model$CholK_1,CholK_K))^2 / c(crossprod(model$CholK_1))
    }
    res$var <- as.vector(var_pred)

  } else {
    # Full covariance matrix
    K_pred <- model$sigma2 * nestedKriging::getCorrMatrix(matrix(X_test, nrow=n_test, ncol=model$d), model$theta, model$kernel) - crossprod(CholK_K)
    if (model$estim_mean==TRUE) {
      K_pred <- K_pred + crossprod(rep(1,length=n_test) - crossprod(model$CholK_1,CholK_K)) / c(crossprod(model$CholK_1))
    }
    res$K <- K_pred
  }

  return(res)
}
