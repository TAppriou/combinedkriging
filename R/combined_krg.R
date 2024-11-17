combined_krg <- function(X_train, Y_train, theta_list, kernel="matern5_2", estim_mean=TRUE, nugget=0, n_tree_levels=NULL, weighting_method="binLOO") {

  model <- list()

  # Initialization of the pre-computed quantities
  model$weighting_method <- weighting_method
  model$estim_mean <- estim_mean
  model$d <- dim(X_train)[2]
  n_train <- dim(X_train)[1]
  model$n_train <- n_train
  model$X_train <- X_train
  model$Y_train <- Y_train

  if (weighting_method=="binLOO") {
    n_comb_max <- 2^(n_tree_levels-1)
  } else {
    n_comb_max <- dim(theta_list)[1]
  }
  model$n_tree_levels <- n_tree_levels
  model$n_comb_max <- n_comb_max
  model$theta_list <- theta_list
  
  model$n_comb_max <- n_comb_max
  model$theta_list <- theta_list
  model$K <- array(0,c(n_train,n_train,n_comb_max))
  model$CholK <- array(0,c(n_train,n_train,n_comb_max))
  model$CholK_1 <- matrix(0,nrow=n_train,ncol=n_comb_max)
  model$mu <- rep(0,length=n_comb_max)
  model$CholK_Y <- matrix(0,nrow=n_train,ncol=n_comb_max)
  model$Kinv_Y <- matrix(0,nrow=n_train,ncol=n_comb_max)
  model$Kinv_diag <- matrix(0,nrow=n_train,ncol=n_comb_max)

  if (kernel=="matern5_2" | kernel=="matern5_2.radial") {
    model$kernel <- "matern5_2.radial"
  } else if (kernel=="exp") {
    model$kernel <- "exp"
  } else if (kernel=="gauss") {
    model$kernel <- "gauss"
  } else if (kernel=="matern3_2" | kernel=="matern3_2.radial") {
    model$kernel <- "matern3_2.radial"
  } else {
    stop("Wrong kernel type")
  }

  # Pre-computed quantities
  for (p in 1:n_comb_max) {
    model$K[,,p] <- nestedKriging::getCorrMatrix(X_train, theta_list[p,], model$kernel) + nugget * diag(n_train)
    model$CholK[,,p] <- chol.default(model$K[,,p])
    model$CholK_1[,p] <- forwardsolve(t(model$CholK[,,p]), matrix(1,nrow=n_train,ncol=1))
    if (estim_mean == TRUE) {
      model$mu[p] <- crossprod(model$CholK_1[,p],forwardsolve(t(model$CholK[,,p]), matrix(Y_train,nrow=n_train,ncol=1))) / crossprod(model$CholK_1[,p])
    }
    model$CholK_Y[,p] <- forwardsolve(t(model$CholK[,,p]), matrix(Y_train-model$mu[p],nrow=n_train,ncol=1))
    model$Kinv_Y[,p] <- backsolve(model$CholK[,,p], model$CholK_Y[,p])
    model$Kinv_diag[,p] <- diag(chol2inv(model$CholK[,,p]))
  }

  if (weighting_method=="MoE") {
    model$beta_list <- weight_MoE(n_comb_max,model)
  } else if (weighting_method=="PoE") {
    model$beta_list <- rep(1,length=n_comb_max)
  } else if (weighting_method=="gPoE") {
    model$beta_list <- weight_gPoE(n_comb_max,model)
  } else if (weighting_method=="LOO") {
    model$beta_list <- weight_LOOCV(n_comb_max,model)
  } else if (weighting_method=="diagLOO") {
    model$beta_list <- weight_LOOCV_diag(n_comb_max,model)

  } else if (weighting_method=="binLOO") {

    # Weights of the sub-models
    beta_nodes <- weight_LOOCV_tree(n_tree_levels, model)
    model$beta_list <- weight_total(n_tree_levels, beta_nodes)

    # Weights of the covariances
    model$alpha_list <- weight_total(n_tree_levels, compute_alpha_tree(beta_nodes, n_tree_levels, model))

    # Global covariance matrix
    sumK <- 0
    for (p in 1:model$n_comb_max) {
      sumK <- sumK + model$alpha_list[p] * model$K[,,p]
    }
    model$sumK <- sumK

    # Cholesky of global covariance matrix
    model$CholsumK <- chol.default(model$sumK)

    # Amplitude of the variance
    model$sigma2 <- pred_s2(model,model$beta_list,model$alpha_list)

  } else {
    stop("Wrong weighting method")
  }

  return(model)
}
