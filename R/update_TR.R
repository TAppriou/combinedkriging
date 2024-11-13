update_TR <- function(d, X_train, Y_train, TR, TR_type="TREGO") {

  if (TR_type=="NoTR") {
    TR$length_TR <- 1
    TR$center_TR <- X_train[which.min(Y_train),]
    TR$upper <- rep(1,d)
    TR$lower <- rep(0,d)
  }


  if (TR_type=="TREGO") {

    # If success
    if (min(Y_train) < TR$bestVal) {

      TR$bestVal <- min(Y_train)
      TR$center_TR <- X_train[which.min(Y_train),]
      TR$length_TR <- min(TR$length_TR * TR$expand_factor, TR$length_max)
      TR$global <- TRUE
      TR$upper <- rep(1,d)
      TR$lower <- rep(0,d)

    } else {

      if (TR$global==FALSE) { #We reduce size only if we fail in local
        TR$length_TR <- TR$length_TR * TR$contract_factor
      }
      if (TR$length_TR < TR$length_min) {
        TR$length_TR <- 0.8
      }

      # If we failed in a global search, we do local search
      if (TR$global==TRUE) {
        TR$global <- FALSE
        for (l in 1:d) {
          TR$upper[l] <- min(TR$center_TR[l] + 0.5*TR$length_TR, 1)
          TR$lower[l] <- max(TR$center_TR[l] - 0.5*TR$length_TR, 0)
        }

        # If we failed in local search, we still go back to global search
      } else {
        TR$global <- TRUE
        TR$upper <- rep(1,d)
        TR$lower <- rep(0,d)
      }

    }

  }

  if (TR_type=="TURBO") {

    # If success
    if (min(Y_train) < TR$bestVal) {

      TR$bestVal <- min(Y_train)
      TR$center_TR <- X_train[which.min(Y_train),]

      TR$count.success <- TR$count.success + 1
      TR$count.failure <- 0
      if (TR$count.success == TR$nsuccess) {
        TR$length_TR <- min(TR$length_TR * TR$expand_factor, TR$length_max)
        TR$count.success <- 0
      }

      # If failure
    } else {
      TR$count.failure <- TR$count.failure + 1
      TR$count.success <- 0
      if (TR$count.failure == TR$nfailure) {
        TR$length_TR <- max(TR$length_TR * TR$contract_factor, TR$length_min)
        if (TR$length_TR < TR$length_min) {
          TR$length_TR <- 0.8
        }
        TR$count.failure <- 0
      }
    }

    # New trust region
    for (l in 1:d) {
      TR$upper[l] <- min(TR$center_TR[l] + 0.5*TR$length_TR, 1)
      TR$lower[l] <- max(TR$center_TR[l] - 0.5*TR$length_TR, 0)
    }

  }

  return(TR)

}
