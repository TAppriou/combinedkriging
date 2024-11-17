\name{combined_krg}
\alias{combined_krg}
\title{Combination of Kriging sub-models with fixed length-scales}
\description{
Build the combination of Kriging sub-models with fixed length-scales.
}
\usage{
combined_krg(X_train, Y_train, theta_list, kernel="matern5_2", estim_mean=TRUE, nugget=0, n_tree_levels, weighting_method="binLOO")
}
\arguments{
  \item{X_train}{The training points. The ith row contains the d input variables corresponding to the ith observation.}
  \item{Y_train}{The function values at the training points. The ith element of the vector is the function value for the ith training point.}
  \item{theta_list}{A matrix containing the length-scales of each sub-model. The ith row contains the d length-scales for the ith sub-model.}
  \item{kernel}{The type of covariance kernel used. The choices are Matérn 5/2 "\code{matern5_2}", Matérn 3/2 "\code{matern3_2}", square exponential "\code{gauss}", or exponential "\code{exp}".}
  \item{theta_init}{An optional vector containing the initial values of the length-scales for the MLE optimization.}
  \item{estim_mean}{Do we use a non-zero constant trend? Default value is \code{TRUE}.}
  \item{nugget}{Optional nugget added to the diagonal of the correlation matrix.}
  \item{n_tree_levels}{For \code{weighting_method=="binLOO"}, number of levels in the binary tree for the combination. The corresponding number of sub-models is \code{2^(n_tree_levels - 1)}.}
  \item{weighting_method}{Weightinh method for the combination of sub-models. The choices are MoE "\code{MoE}", PoE "\code{PoE}", gPoE "\code{gPoE}", LOOCV "\code{LOO}", diagonal LOOCV "\code{diagLOO}", and binary LOOCV "\code{binLOO}" (default method).}
}
\value{
A list of useful pre-computed quantities for prediction containing:
  \item{K}{ An array such that \code{K[,,p]} contains the correltaion matrix for the pth sub-model,}
  \item{CholK}{ An array such that \code{CholK[,,p]} contains the cholesky factor of the correlation matrix for the pth sub-model,}
  \item{mu}{ An array such that \code{mu[p]} contains the constant trend value for the pth sub-model,}
  \item{CholK_1}{ An array such that \code{CholK_1[,p]} contains the inverse Cholesky factor times a vector of ones for the pth sub-model,}
  \item{Kinv_Y}{ An array such that \code{Kinv_Y[,p]} contains the inverse correlation matrix times \code{Y_train - mu} for the pth sub-model,}
  \item{Kinv_diag}{ An array such that \code{Kinv_diag[,p]} contains the diagonal of the inverse correlation matrix for the pth sub-model,}
  \item{beta_list}{ A vector containing the weights of the sub-models in the combination,}
  \item{alpha_list}{ Only available for \code{weighting_method="binLOO"}. A vector containing the weights used to compute the variance of the combination,}
  \item{sumK}{ Only available for \code{weighting_method="binLOO"}. The sum of the correlation matrices of the p sub-models, with weights given in \code{alpha_list},}
  \item{CholsumK}{ Only available for \code{weighting_method="binLOO"}. The Cholesky factor for \code{sumK},}
  \item{sigma2}{ Only available for \code{weighting_method="binLOO"}. Hyperparameter giving the amplitude of the covariance.}
}

\references{
Appriou, T., Rullière, D. and Gaudrie, D., 2023. Combination of optimization-free kriging models for high-dimensional problems. Computational Statistics, pp.1-23.

Appriou, T., Rullière, D. and Gaudrie, D., 2024. High-dimensional Bayesian Optimization with a Combination of Kriging models.
}

\examples{

### An example for the 50D sphere function ###

# Test function: Sphere function
fun_sphere <- function(x) {
  return(sqrt(sum((x-0.5)^2)))
}

d <- 50 #Dimension
n_train <- 200 #Number of training samples
n_test <- 5000 #Number of test points
n_tree_levels <- 5 #Number of levels for the binary tree combinations
n_comb_max <- 2^(n_tree_levels - 1)

set.seed(123) #Set the random seed

# Random training and test points
X_train <- matrix(0,nrow=n_train,ncol=d)
X_test <- matrix(0,nrow=n_test,ncol=d)
for (i in 1:n_train) {
  X_train[i,] <- runif(d, min=0, max=1)
}
for (i in 1:n_test) {
  X_test[i,] <- runif(d, min=0, max=1)
}

# Training and test values
Y_train <- apply(X_train, MARGIN=1, fun_sphere)
Y_test <- apply(X_test, MARGIN=1, fun_sphere)


time1 <- Sys.time()

# Lower and upper bounds for sampling the length-scales
bounds_comb <- bornes_theta(d=d, X=X_train, kernel="matern5_2")

# Sampling the length-scales
theta_list <- sample_theta_entropy(X=X_train, n_comb=n_comb_max, bornes=bounds_comb,
                                   kernel="matern5_2", iso=TRUE)

# Build the combination
model_comb <- combined_krg(X_train=X_train, Y_train=Y_train, theta_list=theta_list,
                           kernel="matern5_2", estim_mean=TRUE, nugget=0,
                           n_tree_levels=n_tree_levels, weighting_method="binLOO")

# Prediction using the combination
prediction <- pred_comb(X_test=X_test, model=model_comb, fullCov=FALSE)

# Compute the R2 coefficient
R2_comb <- 1 - sum((prediction$mean - Y_test)^2)/sum((Y_test-mean(Y_test))^2)

time2 <- Sys.time()
time_comb <- difftime(time2,time1,units="secs")

time3 <- Sys.time()


# Comparison with an ordinary Kriging model using MLE
krg_model <- krg_MLE(X_train=X_train, Y_train=Y_train, bornes_opt=c(0.1,20),
                     kernel="matern5_2", theta_init=NULL, estim_mean=TRUE, nugget=NULL)
prediction_krg <- pred_krg(X_test=X_test, model=krg_model, fullCov=FALSE)

R2_krg <- 1 - sum((prediction_krg$mean - Y_test)^2)/sum((Y_test-mean(Y_test))^2)

time4 <- Sys.time()
time_krg <- difftime(time4,time3,units="secs")

}
