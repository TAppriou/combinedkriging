
# Matern 5/2 covariance function
cov_matern52 <- function(r, theta) {
  dist <- r/theta^2
  return( (1 + sqrt(5*dist) + 5*dist/3) * exp(-sqrt(5*dist)) )
}


# Exponential covariance function
cov_exp <- function(r, theta) {
  dist <- r/theta^2
  return( exp(-sqrt(dist)) )
}


# Square exponential covariance function
cov_gauss <- function(r, theta) {
  dist <- r/theta^2
  return( exp(-0.5*sqrt(dist)) )
}


# Matern 3/2 covariance function
cov_matern32 <- function(r, theta) {
  dist <- r/theta^2
  return( (1 + sqrt(3*dist)) * exp(-sqrt(3*dist)) )
}


# Entropy of the correlation (numerical approximation)
entropy_corr_num <- function(theta, dist_list, kernel="matern5_2") {

  n_samples <- length(dist_list) #Number of observed distances sampled

  #Correlations observed
  if (kernel=="matern5_2") {
    samples_corr <- cov_matern52(dist_list,theta)
  } else if (kernel=="exp") {
    samples_corr <- cov_exp(dist_list,theta)
  } else if (kernel=="gauss") {
    samples_corr <- cov_gauss(dist_list,theta)
  } else if (kernel=="matern3_2") {
    samples_corr <- cov_matern32(dist_list,theta)
  } else {
    stop("Wrong kernel type")
  }

  density_corr <- density(samples_corr, n=1024, from=0, to=1) #Kernel approximation of the density

  # Numerical approximation of the entropy
  entropy <- 0
  for(i in 1:n_samples) {
    ind_min <- which.min(abs(samples_corr[i]-density_corr$x))
    entropy <- entropy - log(density_corr$y[ind_min])
  }
  entropy <- 1/n_samples * entropy

  return(entropy)

}


# Function to sample length-scales using entropy
sample_theta_entropy <- function(X, n_comb, bornes, kernel="matern5_2", iso=TRUE) {

  n_train <- dim(X)[1]
  d <- dim(X)[2]

  # Sampling observed distances
  if (n_train^2 > 10000) {
    n_samples <- 10000  #Number of observed distances sampled
  } else {
    n_samples <- n_train^2
  }
  smp_int <- sample(n_train^2,n_samples)
  pairs_ind <- cbind((smp_int-1) %/% n_train + 1, smp_int %% n_train + 1) #Indices of distances to be sampled
  dist_list <- rep(0, length=n_samples)
  # Sampling
  for (i in 1:n_samples) {
    dist_list[i] <- sum((X[pairs_ind[i,1],]-X[pairs_ind[i,2],])^2)
  }

  theta_list <- matrix(0,nrow=n_comb,ncol=d)

  if (iso==FALSE) {

    for (i in 1:d) {

      theta_values <- seq(bornes[i,1],bornes[i,2],length=200) #Possible values

      # Entropy (numerical approximation)
      h_num <- rep(0,length=200) #Entropy
      for (j in 1:200) {
        h_num[j] <- entropy_corr_num(theta_values[j],dist_list,kernel=kernel)
      }

      # Sampling the length-scales using exp(entropy)
      probas <- exp(h_num)
      theta_list[,i] <- sample(x = theta_values, size = n_comb,replace = TRUE, prob=probas)
    }

  } else {

    theta_values <- seq(min(bornes[,1]),max(bornes[,2]),length=200) #Possible values

    # Entropy (numerical approximation)
    h_num <- rep(0,length=200) #Entropy
    for (j in 1:200) {
      h_num[j] <- entropy_corr_num(theta_values[j],dist_list,kernel=kernel)
    }

    # Sampling the length-scales using exp(entropy)
    probas <- exp(h_num)

    for (i in 1:d) {
      theta_list[,i] <- sample(x = theta_values, size = n_comb,replace = TRUE, prob=probas)
    }

  }

  return(theta_list)

}
