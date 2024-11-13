\name{optim_EI}
\alias{optim_EI}
\title{Optimize the Expected Improvement}
\description{
Optimize the Expected Improvement for the combination of for ordinary Kriging.
}
\usage{
optim_EI(model, upper, lower, full.exploitation=FALSE, combination=TRUE, optimizer="cmaes")
}
\arguments{
  \item{model}{A combination of Kriging sub-models obtained with the function \code{combined_krg}, or an ordinary Kriging model obtained with \code{krg_MLE}.}
  \item{upper}{A vector giving the upper bound in each dimension of the box where the EI is maximized.}
  \item{lower}{A vector giving the lower bound in each dimension of the box where the EI is maximized.}
  \item{full.exploitation}{\code{FALSE} if we maximize the EI (default). \code{TRUE} if we minimize the Kriging mean instead.}
  \item{combination}{\code{TRUE} for a combined Kriging model, \code{FALSE} for an ordinary Kriging model.}
  \item{optimizer}{The optimizer for the EI. Only "cmaes" is implemented in this version.}
}
\value{
A list containing:
\item{par}{The point which maximizes the EI,}
\item{value}{-EI at this location.}
}
