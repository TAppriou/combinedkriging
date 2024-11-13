\name{bornes_theta}
\alias{bornes_theta}
\title{Bounds for the length-scale sampling}
\description{
Compute the bounds in which the length-scales will be sampled.
}
\usage{
bornes_theta(d, X, kernel="matern5_2")
}
\arguments{
  \item{d}{The dimension of the problem.}
  \item{X}{The design of experiments. The ith row contains of the d input variables corresponding to the ith observation.}
  \item{kernel}{The type of covariance kernel used. The choices are Matérn 5/2 "\code{matern5_2}", Matérn 3/2 "\code{matern3_2}", square exponential "\code{gauss}", or exponential "\code{exp}".}
}
\value{The lower and upper bounds for the length-scale sampling.}


