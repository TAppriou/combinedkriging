\name{weight_LOOCV}
\alias{weight_LOOCV}
\title{Compute the LOOCV weights}
\description{
Compute the sub-models weights using the Leave-one-out Cross-Validation (LOOCV) method.
}
\usage{
weight_LOOCV(n_comb, model, nugget=0)
}
\arguments{
  \item{n_comb}{Number of sub-models in the combination}
  \item{model}{A combination of Kriging sub-models obtained with the function \code{combined_krg}.}
  \item{nugget}{Optional nugget added in the diagonal of the LOOCV residuals matrix to improve numerical stability.}
}
\value{A vector containing the d weights for the sub-models.}

\references{
Appriou, T., Rullière, D. and Gaudrie, D., 2023. Combination of optimization-free kriging models for high-dimensional problems. Computational Statistics, pp.1-23.
}
