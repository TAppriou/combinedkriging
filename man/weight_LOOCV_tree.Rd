\name{weight_LOOCV_tree}
\alias{weight_LOOCV_tree}
\title{Compute the LOOCV weights on a binary tree}
\description{
Compute the weights at each node of a binary tree using Leave-one-out Cross-Validation (LOOCV).
}
\usage{
weight_LOOCV_tree(n_tree_levels, model)
}
\arguments{
  \item{n_tree_levels}{Number of levels in the binary tree. The number of nodes is \code{2^n_tree_levels - 1}, and the corresponding number of sub-models is \code{2^(n_tree_levels - 1)}.}
  \item{model}{A combination of Kriging sub-models obtained with the function \code{combined_krg}.}
}
\value{A vector containing the weights for each node of the binary tree.}

\references{
Appriou, T., Rullière, D. and Gaudrie, D., 2023. Combination of optimization-free kriging models for high-dimensional problems. Computational Statistics, pp.1-23.

Appriou, T., Rullière, D. and Gaudrie, D., 2024. High-dimensional Bayesian Optimization with a Combination of Kriging models.
}
