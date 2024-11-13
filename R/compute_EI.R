compute_EI <- function(X, model, full.exploitation=FALSE, combination=TRUE) {

  # Ordinary Kriging (no combination)
  if (combination==FALSE) {

    prediction <- pred_krg(X, model)
    # Prediction mean and variance
    krg_mean <- prediction$mean
    if (sum((prediction$var < 0)) > 0) {
      stop("Error : negative variance !")
    } else {
      krg_sd <- sqrt(prediction$var)
    }

    # Combination of Kriging sub-models
  } else {
    prediction <- pred_comb(X, model)
    # Prediction mean and variance
    krg_mean <- prediction$mean
    krg_sd <- sqrt(prediction$var)
    if (sum((prediction$var < 0)) > 0) {
      stop("Error : negative variance")
    } else {
      krg_sd <- sqrt(prediction$var)
    }
  }

  # Best current value
  current_best <- min(model$Y_train)

  # EI
  improv <- (current_best - krg_mean) / krg_sd
  improv.prob <- pnorm(improv)
  improv.dens <- dnorm(improv)
  if (!full.exploitation) {
    res <- (current_best - krg_mean) * improv.prob + krg_sd * improv.dens
  } else {
    res <- (current_best - krg_mean)
  }

  return(-res)

}
