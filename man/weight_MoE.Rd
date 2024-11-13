\name{weight_MoE}
\alias{weight_MoE}
\title{Compute the MoE weights}
\description{
Compute the sub-models weights using the Mixture of Experts (MoE) method.
}
\usage{
weight_MoE(n_comb, model)
}
\arguments{
  \item{n_comb}{Number of sub-models in the combination}
  \item{model}{A combination of Kriging sub-models obtained with the function \code{combined_krg}.}
}
\value{A vector containing the d weights for the sub-models.}

\references{
Appriou, T., Rullière, D. and Gaudrie, D., 2023. Combination of optimization-free kriging models for high-dimensional problems. Computational Statistics, pp.1-23.
}
