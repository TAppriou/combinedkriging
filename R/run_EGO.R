run_EGO <- function(n_steps, fun, ..., X_train, Y_train, model, bornes_theta_opt=c(0.1,20), estim_mean=TRUE, combination=TRUE, n_tree_levels=NULL, optimizer="cmaes", TR_type="TREGO", nugget=1e-6, log_file=NULL) {

  if (optimizer!="cmaes") {
    stop("Only cmaes implemented for this version")
  }

  # Initialize trust region
  TR <- initialize_TR(model$d, X_train, Y_train, TR_type)

  # History of the trust region
  history_center_TR <- matrix(0, nrow=n_steps+1, ncol=model$d)
  history_length_TR <- rep(0, length=n_steps+1)
  history_center_TR[1,] <- TR$center_TR
  history_length_TR[1] <- TR$length_TR

  # Ordinary Kriging (no combination)
  if (combination==FALSE) {

    time_MLE <- 0
    time_EI <- 0

    # EGO loop
    i <- 1
    while (i < n_steps) {

      full.exploitation <- FALSE

      # Number of iter before updating the trust region
      if (TR$global == TRUE) {
        n_iter <- TR$n_global
      } else {
        n_iter <- TR$n_local
      }

      for (j in 1:n_iter) {

        if (i < n_steps) {

          if (!is.null(log_file)) {
            msg <- paste("(Ordinary Kriging) Starting EGO loop",i)
            file_name <- paste(log_file,"/OK_iteration",i,".txt",sep="")
            cat(msg, file=file_name)
          } else {
            message(paste("(Ordinary Kriging) Starting EGO loop",i))
          }

          time1 <- Sys.time()

          # New point with EI
          max_EI <- optim_EI(model, TR$upper, TR$lower, full.exploitation, combination=FALSE, optimizer)

          # Add the new point to the DoE
          X_next <- max_EI$par
          X_train <- rbind(X_train, X_next, deparse.level=0)
          newVal <- fun(X_next, ...)
          Y_train <- c(Y_train, newVal)

          time2 <- Sys.time()
          time_EI <- time_EI + difftime(time2,time1,units="mins")

          # Re-optimize the hyperparameters
          model <- krg_MLE(X_train, Y_train, bornes_theta_opt, model$kernel, theta_init=model$theta, estim_mean, nugget)

          time3 <- Sys.time()
          time_MLE <- time_MLE + difftime(time3,time2,units="mins")

          i <- i+1

        }

      }

      # Update the trust region
      TR <- update_TR(model$d, X_train, Y_train, TR, TR_type)
      history_center_TR[i,] <- TR$center_TR
      history_length_TR[i] <- TR$length_TR

    }

    # One full exploitation at the end
    full.exploitation <- TRUE
    max_EI <- optim_EI(model, rep(1,model$d), rep(0,model$d), full.exploitation, combination=FALSE, optimizer)
    X_next <- max_EI$par
    X_train <- rbind(X_train, X_next, deparse.level=0)
    newVal <- fun(X_next, ...)
    Y_train <- c(Y_train, newVal)
    model <- krg_MLE(X_train, Y_train, bornes_theta_opt, model$kernel, theta_init=model$theta, estim_mean, nugget)

    time_list <- list(time_EI=time_EI, time_MLE=time_MLE)


    # Combination of Kriging sub-models
  } else {

    if (model$weighting_method!="binLOO") {
      stop("EGO only implemented for the binLOO weighting method")
    }

    time_EI <- 0
    time_combination <- 0

    # EGO loop
    i <- 1
    while (i < n_steps) {

      full.exploitation <- FALSE

      # Number of iter before updating the trust region
      if (TR$global == TRUE) {
        n_iter <- TR$n_global
      } else {
        n_iter <- TR$n_local
      }

      for (j in 1:n_iter) {

        if (i < n_steps) {

          if (!is.null(log_file)) {
            msg <- paste("(Combination) Starting EGO loop",i)
            file_name <- paste(log_file,"/comb_iteration",i,".txt",sep="")
            cat(msg, file=file_name)
          } else {
            message(paste("(Combination) Starting EGO loop",i))
          }

          time1 <- Sys.time()

          # New point with EI
          max_EI <- optim_EI(model, TR$upper, TR$lower, full.exploitation, combination=TRUE, optimizer)

          # Add the new point to the DoE
          X_next <- max_EI$par
          X_train <- rbind(X_train, X_next, deparse.level=0)
          newVal <- fun(X_next, ...)
          Y_train <- c(Y_train, newVal)

          time2 <- Sys.time()
          time_EI <- time_EI + difftime(time2,time1,units="mins")

          # Re-build the sub-models
          model <- combined_krg(X_train, Y_train, model$theta_list, model$kernel, estim_mean, nugget, model$n_tree_levels, model$weighting_method)

          time3 <- Sys.time()
          time_combination <- time_combination + difftime(time3,time2,units="mins")

          i <- i+1

        }

      }

      # Update the trust region
      TR <- update_TR(model$d, X_train, Y_train, TR, TR_type)
      history_center_TR[i,] <- TR$center_TR
      history_length_TR[i] <- TR$length_TR

    }

    # One full exploitation at the end
    full.exploitation <- TRUE
    max_EI <- optim_EI(model, TR$upper, TR$lower, full.exploitation, combination=TRUE, optimizer)
    X_next <- max_EI$par
    X_train <- rbind(X_train, X_next, deparse.level=0)
    newVal <- fun(X_next, ...)
    Y_train <- c(Y_train, newVal)
    model <- combined_krg(X_train, Y_train, model$theta_list, model$kernel, estim_mean, nugget, model$n_tree_levels, model$weighting_method)

    time_list <- list(time_EI=time_EI, .time_combination=time_combination)

  }

  return( list(model=model, X=X_train, Y=Y_train, time=time_list, history_center_TR=history_center_TR, history_length_TR=history_length_TR) )

}
