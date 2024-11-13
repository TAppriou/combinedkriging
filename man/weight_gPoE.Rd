\name{weight_gPoE}
\alias{weight_gPoE}
\title{Compute the gPoE weights}
\description{
Compute the inner weights for the generalized Product of Experts (gPoE) weighting method.
}
\usage{
weight_gPoE(n_comb, model)
}
\arguments{
  \item{n_comb}{Number of sub-models in the combination}
  \item{model}{A combination of Kriging sub-models obtained with the function \code{combined_krg}.}
}
\value{A vector containing the d inner gPoE weights.}

\references{
Appriou, T., Rullière, D. and Gaudrie, D., 2023. Combination of optimization-free kriging models for high-dimensional problems. Computational Statistics, pp.1-23.
}
