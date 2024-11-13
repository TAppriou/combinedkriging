sample_theta_copula <- function(d, n_comb, bornes, param_copula=2) {

  theta_list <- matrix(0,nrow=n_comb,ncol=d)
  mycopula <- copula::gumbelCopula(param=param_copula, dim=d)
  smp_copula <- copula::rCopula(n_comb, mycopula)
  for (i in 1:d) {
    theta_list[,i] <- bornes[i,2] - smp_copula[,i] * (bornes[i,2]-bornes[i,1])
  }

  return(theta_list)

}
