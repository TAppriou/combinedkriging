\name{compute_alpha_tree}
\alias{compute_alpha_tree}
\title{Compute the weights for the variance of the combination on a binary tree}
\description{
Compute the weights at each node of a binary tree for the variance of the combination, based on the weights of the combination on the tree.
}
\usage{
compute_alpha_tree(w_nodes, n_tree_levels, model)
}
\arguments{
  \item{w_nodes}{The weights of the combination on each node of the binary tree,}
  \item{n_tree_levels}{Number of levels in the binary tree. The number of nodes is \code{2^n_tree_levels - 1}, and the corresponding number of sub-models is \code{2^(n_tree_levels - 1)}.}
  \item{model}{A combination of Kriging sub-models obtained with the function \code{combined_krg}.}
}
\value{A vector containing the variance weights for each node of the binary tree.}

\references{

Appriou, T., Rullière, D. and Gaudrie, D., 2024. High-dimensional Bayesian Optimization with a Combination of Kriging models.
}
