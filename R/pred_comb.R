pred_comb <- function(X_test, model, fullCov=FALSE) {

  if (model$weighting_method!="binLOO") {
    stop("This function is only for weighting_method=='binLOO'")
  }

  res <- list()

  if (is.matrix(X_test)) {
    n_test <- dim(X_test)[1]
  } else {
    n_test <- 1
  }

  n_train <- model$n_train
  n_comb_max <- model$n_comb_max

  sum_kXx <- 0
  # Mean Prediction for each submodel
  Y_pred_submodels <- matrix(0,nrow=model$n_comb_max, ncol=n_test) #Initialization
  for (p in 1:n_comb_max) {
    K_X_Xtest <- nestedKriging::getCrossCorrMatrix(model$X_train, matrix(X_test, nrow=n_test, ncol=model$d), model$theta_list[p,], model$kernel)
    CholK_K <- forwardsolve(t(model$CholK[,,p]), K_X_Xtest)
    Y_pred_submodels[p,] <- matrix(model$mu[p], nrow=n_test, ncol=1) + crossprod(CholK_K, model$CholK_Y[,p])
    sum_kXx <- sum_kXx + model$alpha_list[p] * K_X_Xtest
  }

  # Mean Prediction for the combination
  Y_pred <- colSums(model$beta_list * Y_pred_submodels)
  res$mean <- as.vector(Y_pred)

  CholsumK_sumk <- forwardsolve(t(model$CholsumK),sum_kXx)

  if (fullCov == FALSE) {
    var_pred <- model$sigma2 * (rep(sum(model$alpha_list),length=n_test) - colSums(CholsumK_sumk^2))
    res$var <- as.vector(var_pred)

    # Full covariance matrix
  } else {
    K_pred <- matrix(0, nrow=n_test, ncol=n_test)
    for (p in 1:n_comb_max) {
      K_pred <- K_pred + model$alpha_list[p] * nestedKriging::getCorrMatrix(matrix(X_test, nrow=n_test, ncol=model$d), model$theta_list[p,], model$kernel)
    }
    K_pred <- model$sigma2 * (K_pred - crossprod(CholsumK_sumk))
    res$K <- K_pred
  }


  return(res)
}
