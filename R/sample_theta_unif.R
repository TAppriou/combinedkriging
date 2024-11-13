sample_theta_unif <- function(d, n_comb, bornes) {

  theta_list <- matrix(0,nrow=n_comb,ncol=d)
  for (i in 1:d) {
    theta_list[,i] <- runif(n_comb, min=bornes[i,1], max=bornes[i,2])
  }

  return(theta_list)

}
