bornes_theta <- function(d, X, kernel="matern5_2") {

  bornes <- matrix(0,nrow=d,ncol=2)

  # Typical distances in the DoE
  kappa <- 9/5
  if (d>6) {
    r_min <- sqrt(2*d - 1.96*sqrt(2*(kappa+1)*d))
  } else {
    r_min <- 1
  }
  r_max <- sqrt(2*d + 1.96*sqrt(2*(kappa+1)*d))

  # Kernel impact
  if (kernel=="matern5_2") {
    theta_m <- 0.22991
    theta_p <- 2.4411
  } else if (kernel=="exp") {
    theta_m <- 0.14861
    theta_p <- 3.7732
  } else if (kernel=="gauss") {
    theta_m <- 0.29255
    theta_p <- 1.9641
  } else if (kernel=="matern3_2") {
    theta_m <- 0.20652
    theta_p <- 2.7384
  } else {
    stop("Wrong kernel type")
  }

  for (i in 1:d) {
    s <- sd(X[,i])
    bornes[i,1] <- s*r_min*theta_m
    bornes[i,2] <- s*r_max*theta_p
  }

  return(bornes)
}
