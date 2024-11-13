weight_LOOCV_diag <- function(n_comb, model) {

  diagA <- rep(0,length=n_comb)
  for (i in 1:model$n_train) {
    diagA <- diagA + (model$Kinv_Y[i,1:n_comb]/model$Kinv_diag[i,1:n_comb])^2
  }

  return( c((1/diagA)/sum(1/diagA)) )

}
