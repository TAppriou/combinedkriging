pred_comb_mean <- function(X_test, model) {

  res <- list()

  if (is.matrix(X_test)) {
    n_test <- dim(X_test)[1]
  } else {
    n_test <- 1
  }

  n_train <- model$n_train
  n_comb_max <- model$n_comb_max

  # Mean Prediction for each submodel
  Y_pred_submodels <- matrix(0,nrow=model$n_comb_max, ncol=n_test) #Initialization
  var_submodels <- matrix(0,nrow=model$n_comb_max, ncol=n_test)
  for (p in 1:n_comb_max) {
    K_X_Xtest <- nestedKriging::getCrossCorrMatrix(model$X_train, matrix(X_test, nrow=n_test, ncol=model$d), model$theta_list[p,], model$kernel)
    CholK_K <- forwardsolve(t(model$CholK[,,p]), K_X_Xtest)
    Y_pred_submodels[p,] <- matrix(model$mu[p], nrow=n_test, ncol=1) + crossprod(CholK_K, model$CholK_Y[,p])
    var_submodels[p,] <- rep(1,length=n_test) - colSums(CholK_K^2)
  }

  # Mean Prediction for the combination
  if (model$weighting_method=="PoE" | model$weighting_method=="gPoE") {
    Y_pred <- colSums(model$beta_list/var_submodels * Y_pred_submodels) / (model$beta_list %*% (1/var_submodels))
  } else {
    Y_pred <- colSums(model$beta_list * Y_pred_submodels)
  }

  res$mean <- as.vector(Y_pred)

  return(res)
}
