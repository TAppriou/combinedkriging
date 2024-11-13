initialize_TR <- function(d, X_train, Y_train, TR_type="TREGO") {

  # Initialize list
  TR <- list()

  # General
  TR$expand_factor <- 2
  TR$contract_factor <- 0.5
  TR$length_TR <- 0.8
  TR$length_max <- Inf
  TR$length_min <- 2^-7
  TR$bestVal <- min(Y_train)
  TR$center_TR <- X_train[which.min(Y_train),]
  TR$upper <- rep(1,d)
  TR$lower <- rep(0,d)

  if (TR_type=="NoTR") {
    TR$global <- TRUE
    TR$n_global <- 1
  }

  if (TR_type=="TREGO") {
    TR$global <- TRUE
    TR$n_global <- 1
    TR$n_local <- 4
  }

  if (TR_type=="TURBO") {
    TR$global <- FALSE
    TR$n_local <- 1
    TR$nsuccess <- 3
    TR$nfailure <- min(d,10)
    TR$count.success <- 0
    TR$count.failure <- 0
    for (l in 1:d) {
      TR$upper[l] <- min(TR$center_TR[l] + 0.5*TR$length_TR, 1)
      TR$lower[l] <- max(TR$center_TR[l] - 0.5*TR$length_TR, 0)
    }
  }

  return(TR)

}
