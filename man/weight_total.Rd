\name{weight_total}
\alias{weight_total}
\title{Compute the sub-models weights for a binary tree}
\description{
Compute the sub-models weights from the weights on each nodes for a binary tree structure.
}
\usage{
weight_total(n_tree_levels, w)
}
\arguments{
  \item{n_tree_levels}{Number of levels in the binary tree. The number of nodes is \code{2^n_tree_levels - 1}, and the corresponding number of sub-models is \code{2^(n_tree_levels - 1)}.}
  \item{w}{A vector containing the weights for each node of the binary tree.}
}
\value{A vector containing the weights of each sub-model.}

\references{
Appriou, T., Rullière, D. and Gaudrie, D., 2023. Combination of optimization-free kriging models for high-dimensional problems. Computational Statistics, pp.1-23.

Appriou, T., Rullière, D. and Gaudrie, D., 2024. High-dimensional Bayesian Optimization with a Combination of Kriging models.
}
