\name{sample_theta_copula}
\alias{sample_theta_copula}
\title{Sampling of the length-scales with copulas}
\description{
Sample the length-scales using Gumbel copulas with uniform marginals.
}
\usage{
sample_theta_copula(d, n_comb, bornes, param_copula=2)
}
\arguments{
  \item{d}{The dimension of the problem.}
  \item{n_comb}{Number of length-scale vectors to be sampled.}
  \item{bornes}{The lower and upper bounds for the sampling.}
  \item{param_copula}{The parameter of the Gumbel copula. Must be greater or equal to 1.}
}
\value{The length-scales samples. The ith row contains the d length-scales for the ith sample.}


