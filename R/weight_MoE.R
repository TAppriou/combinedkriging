weight_MoE <- function(n_comb,model) {

  L <- rep(0,length=n_comb)
  for (p in 1:n_comb) {
    s2 <- crossprod(model$CholK_Y[,p])/model$n_train
    term1 <- -model$n_train/2 * log(2*pi)
    term2 <- -model$n_train/2*log(s2)
    term3 <- -1/2*2*sum(log(diag(model$CholK[,,p])))
    term4 <- -1/2*model$n_train
    L[p] <- term1+term2+term3+term4
  }
  mL <- max(L)

  return( exp(L-mL)/sum(exp(L-mL)) )

}
