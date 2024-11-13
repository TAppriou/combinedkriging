weight_LOOCV <- function(n_comb, model, nugget=0) {

  A <- matrix(0,nrow=n_comb,ncol=n_comb) + nugget*diag(n_comb)
  for (i in 1:model$n_train) {
    A <- A + tcrossprod(model$Kinv_Y[i,1:n_comb]/model$Kinv_diag[i,1:n_comb])
  }
  Ainv <- solve(A, matrix(1,nrow=n_comb,ncol=1))

  return( c(Ainv/c(matrix(1,nrow=1,ncol=n_comb) %*% Ainv)) )

}
