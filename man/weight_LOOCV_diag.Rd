\name{weight_LOOCV_diag}
\alias{weight_LOOCV_diag}
\title{Compute the diagonal LOOCV weights}
\description{
Compute the sub-models weights using the diagonal Leave-one-out Cross-Validation (LOOCV) method.
}
\usage{
weight_LOOCV_diag(n_comb, model)
}
\arguments{
  \item{n_comb}{Number of sub-models in the combination}
  \item{model}{A combination of Kriging sub-models obtained with the function \code{combined_krg}.}
}
\value{A vector containing the d weights for the sub-models.}

\references{
Appriou, T., Rullière, D. and Gaudrie, D., 2023. Combination of optimization-free kriging models for high-dimensional problems. Computational Statistics, pp.1-23.
}
