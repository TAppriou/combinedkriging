# Transformation of the weight to the simplex
w_simplex <- function(alpha) {
  w <- rep(0,length=length(alpha)+1)
  alpha_ex <- c(alpha,1)
  w[1] <- alpha_ex[1]
  for (i in 2:length(w)) {
    w[i] <- alpha_ex[i]
    for (k in 1:(i-1)) {
      w[i] <- w[i]*(1-alpha_ex[k])
    }
  }
  return(w)
}


# Gradient of the simplex transformation
dw_simplex <- function(alpha,j) {
  dw <- rep(0,length=length(alpha)+1) #Grad of the transformation
  alpha_ex <- c(alpha,1)
  for (i in 1:length(dw)) {
    if (j>i) {
      dw[i] <- 0
    }
    if (j<i) {
      dw[i] <- -alpha_ex[i]
      for (k in (1:(i-1))[-j]) {
        dw[i] <- dw[i]*(1-alpha_ex[k])
      }
    }
    if (j==i) {
      dw[i] <- 1
      if (i!=1) {
        for (k in 1:(i-1)) {
          dw[i] <- dw[i]*(1-alpha_ex[k])
        }
      }
    }
  }
  return(dw)
}



# Weights of the gPoE method for the combination (optimized with the LOOCV error)
weight_gPoE <- function(n_comb, model) {

  # LOOCV error
  LOOCV_error <- function(alpha) {
    w <- w_simplex(alpha) #Transformation to the simplex
    return( sum(( (model$Kinv_Y[,1:n_comb] %*% w)/(model$Kinv_diag[,1:n_comb] %*% w) )^2) ) #LOOCV error
  }

  # Gradient of the LOOCV error
  grad_LOOCV_error <- function(alpha) {
    grad <- rep(0,length=n_comb-1)
    w <- w_simplex(alpha)
    for (j in 1:n_comb-1) {
      dw <- dw_simplex(alpha,j) #Grad of the transformation
      grad[j] <- sum( (model$Kinv_Y[,1:n_comb] %*% w)/((model$Kinv_diag[,1:n_comb] %*% w)^2) * ( (model$Kinv_Y[,1:n_comb] %*% dw) - (model$Kinv_Y[,1:n_comb] %*% w)*(model$Kinv_diag[,1:n_comb] %*% dw)/(model$Kinv_diag[,1:n_comb] %*% w) ) )
    }
    return(grad)
  }

  # Optimization
  min_err <- 1e99
  for (i in 1:1) {
    alpha_init <- matrix( runif(20*(n_comb-1), min=0, max=1), nrow=20, ncol=n_comb-1 ) # 20 Initial random weights for the optim
    alpha_init <- alpha_init[ which.min(apply(X=alpha_init, MARGIN=1, FUN=LOOCV_error)), ]
    opt <- optim(par=alpha_init, LOOCV_error, grad_LOOCV_error, method="L-BFGS-B",
                 lower=rep(0,length=n_comb-1), upper=rep(1,length=n_comb-1) )
    if (opt$value < min_err) {
      min_err <- opt$value
      res <- opt
    }
  }

  # Weights
  return( w_simplex(res$par) )
}
