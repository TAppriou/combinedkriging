\name{krg_MLE}
\alias{krg_MLE}
\title{Function to build an  ordinary Kriging model}
\description{
Build an ordinary Kriging model using maximum likelihood estimation and save usuful pre_computed quantities. Based on the package DiceKriging.
}
\usage{
krg_MLE(X_train, Y_train, bornes_opt, kernel="matern5_2", theta_init=NULL, estim_mean=TRUE, nugget=1e-6)
}
\arguments{
  \item{X_train}{The training points. The ith row contains the d input variables corresponding to the ith observation.}
  \item{Y_train}{The function values at the training points. The ith element of the vector is the function value for the ith training point.}
  \item{bornes_opt}{Bounds used for the length-scale MLE optimization.}
  \item{kernel}{The type of covariance kernel used. The choices are Matérn 5/2 "\code{matern5_2}", Matérn 3/2 "\code{matern3_2}", square exponential "\code{gauss}", or exponential "\code{exp}".}
  \item{theta_init}{An optional vector containing the initial values of the length-scales for the MLE optimization.}
  \item{estim_mean}{Do we use a non-zero constant trend? Default value is \code{TRUE}}.
  \item{nugget}{Nugget added to the diagonal of the correlation matrix. By default, \code{nugget=1e-6} to avoid numerical issues while inverting the matrix.}
}
\value{
A list of useful quantities for prediction containing:
  \item{CholK}{ Cholesky factor of the covariance matrix,}
  \item{mu}{ Constant trend value,}
  \item{CholK_1}{ Inverse Cholesky factor times a vector of ones,}
  \item{CholK_Y}{ Inverse Cholesky factor times \code{Y_train - mu},}
  \item{sigma2}{Hyperparameter giving the amplitude of the covariance,}
  \item{Covstruct}{ Object from the package DiceKriging containing the covariance structure.}
}
\references{
Roustant, O., Ginsbourger, D. and Deville, Y., 2012. \emph{DiceKriging, DiceOptim: Two R packages for the analysis of computer experiments by kriging-based metamodeling and optimization}. Journal of statistical software, 51, pp.1-55.
}
