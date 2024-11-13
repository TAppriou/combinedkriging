pred_s2 <- function(model, beta_list, alpha_list) {

  n_train <- model$n_train

  # Inverse of the global covariance matrix
  Ktotinv_diag <- diag(chol2inv(model$CholsumK))

  # LOO residuals of the combination
  eLOO <- rep(0,length=n_train)
  for (i in 1:n_train) {
    tmp <- 0
    for (p in 1:model$n_comb_max) {
      tmp <- tmp + beta_list[p] * model$Kinv_Y[i,p]/model$Kinv_diag[i,p] #LOO of the combination
    }
    eLOO[i] <- tmp*sqrt(Ktotinv_diag[i]) #Normalization by the global GP variance
  }

  # Amplitude of the variance hyperparameter
  model$sigma2 <- ( (quantile(eLOO,0.75)-quantile(eLOO,0.25))/(2*qnorm(0.75,mean=0,sd=1)) )^2

}
