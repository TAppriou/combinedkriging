\name{compute_EI}
\alias{compute_EI}
\title{Compute the Expected Improvement}
\description{
Compute the Expected Improvement for the combination of for ordinary Kriging.
}
\usage{
compute_EI(X, model, full.exploitation=FALSE, combination=TRUE)
}
\arguments{
  \item{X}{The locations where we compute the EI.}
  \item{model}{A combination of Kriging sub-models obtained with the function \code{combined_krg}, or an ordinary Kriging model obtained with \code{krg_MLE}.}
  \item{full.exploitation}{\code{FALSE} if we want the EI (default). \code{TRUE} gives the Kriging mean instead.}
  \item{combination}{\code{TRUE} for the EI of a combined Kriging model, \code{FALSE} for the EI of an ordinary Kriging model.}
}
\value{
A vector giving the value of -EI (to be minimized) at each location in \code{X}.
}
