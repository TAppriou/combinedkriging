\name{pred_krg}
\alias{pred_krg}
\title{Prediction using an ordinary Kriging model}
\description{
Compute the mean and variance predictions for an ordinary Kriging model.
}
\usage{
pred_krg(X_test, model, fullCov=FALSE)
}
\arguments{
  \item{X_test}{The locations where we compute the prediction. The ith row contains the d input variables corresponding to the ith location.}
  \item{model}{A kriging model built with the fonction \code{krg_MLE}.}
  \item{fullCov}{\code{FALSE} if we do not need the full covariance matrix between test points (default).}
}
\value{
A list containing the mean and variance prediction:
  \item{mean}{ A vector giving the Kriging mean at each test location.}
  \item{var}{ A vector giving the Kriging variance at each test location.}
  \item{K}{ If \code{fullCov=TRUE}, the covariance matrix between the test points.}
}

