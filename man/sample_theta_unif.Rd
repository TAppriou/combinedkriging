\name{sample_theta_unif}
\alias{sample_theta_unif}
\title{Uniform sampling of the length-scales}
\description{
Sample the length-scales uniformly i.i.d between the bounds.
}
\usage{
sample_theta_unif(d, n_comb, bornes)
}
\arguments{
  \item{d}{The dimension of the problem.}
  \item{n_comb}{Number of length-scale vectors to be sampled.}
  \item{bornes}{The lower and upper bounds for the sampling.}
}
\value{The length-scales samples. The ith row contains the d length-scales for the ith sample.}

