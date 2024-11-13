\name{sample_theta_entropy}
\alias{sample_theta_entropy}
\title{Sampling of the length-scales with the entropy of the correlations}
\description{
Sample the length-scales between the bounds using the entropy of the correlations.
}
\usage{
sample_theta_entropy(X, n_comb, bornes, kernel="matern5_2", iso=TRUE)
}
\arguments{
  \item{X}{The design of experiments. The ith row contains of the d input variables corresponding to the ith observation.}
  \item{n_comb}{Number of length-scale vectors to be sampled.}
  \item{bornes}{The lower and upper bounds for the sampling.}
  \item{kernel}{The type of covariance kernel used. The choices are Matérn 5/2 "\code{matern5_2}", Matérn 3/2 "\code{matern3_2}", square exponential "\code{gauss}", or exponential "\code{exp}".}
  \item{iso}{Do we use the same entropy approximation to sample in every dimension? It saves time when this is the case.}
}
\value{The length-scales samples. The ith row contains the d length-scales for the ith sample.}

