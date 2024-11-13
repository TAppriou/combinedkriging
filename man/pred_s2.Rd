\name{pred_s2}
\alias{pred_s2}
\title{Compute the hyperparameter giving the amplitude of the covariance}
\description{
Compute the hyperparameter giving the amplitude of the covariance based on the weights of the sub-models and of the variance. The amplitude is estimated by fitting the emprical inter-quartile distance of the normalized LOOCV residuals to that of a standard normal distribution.
}
\usage{
pred_s2(model, beta_list, alpha_list)
}
\arguments{
  \item{model}{A combination of Kriging sub-models obtained with the function \code{combined_krg},}
  \item{beta_list}{A vector containing the weights of the sub-models in the combination}
  \item{alpha_list}{A vector containing the weights used to compute the variance of the combination}
}
\value{The amplitude of the covariance.}

\references{

Appriou, T., Rullière, D. and Gaudrie, D., 2024. High-dimensional Bayesian Optimization with a Combination of Kriging models.
}
