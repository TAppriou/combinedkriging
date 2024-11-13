\name{update_TR}
\alias{update_TR}
\title{Update the Trust Region}
\description{
Update the Trust Region.
}
\usage{
update_TR(d, X_train, Y_train, TR, TR_type="TREGO")
}
\arguments{
  \item{d}{The dimension of the problem.}
  \item{X_train}{The training points. The ith row contains the d input variables corresponding to the ith observation.}
  \item{Y_train}{The function values at the training points. The ith element of the vector is the function value for the ith training point.}
  \item{TR}{A trust region created with the function \code{initialize_TR}.}
  \item{TR_type}{The type of Trust region used. Options are no trust region "\code{NoTR}", TREGO implementation of trust regions "\code{TREGO}", TURBO implementation of trust regions "\code{TURBO}".}
}
\value{
The updated Trust Region list.
}

\references{

Diouane, Y., Picheny, V., Riche, R.L. and Perrotolo, A.S.D., 2023. TREGO: a trust-region framework for efficient global optimization. Journal of Global Optimization, 86(1), pp.1-23.

Eriksson, D., Pearce, M., Gardner, J., Turner, R.D. and Poloczek, M., 2019. Scalable global optimization via local Bayesian optimization. Advances in neural information processing systems, 32.

}
