\name{pred_comb}
\alias{pred_comb}
\title{Prediction using the combined Kriging model}
\description{
Compute the mean and variance predictions for the combination of Kriging sub-models.
}
\usage{
pred_comb(X_test, model, fullCov=FALSE)
}
\arguments{
  \item{X_test}{The locations where we compute the prediction. The ith row contains the d input variables corresponding to the ith location.}
  \item{model}{A combination of Kriging sub-models obtained with the function \code{combined_krg}.}
  \item{fullCov}{\code{FALSE} if we do not need the full covariance matrix between test points (default).}
}
\value{
A list containing the mean and variance prediction:
  \item{mean}{ A vector giving the Kriging mean at each test location.}
  \item{var}{ A vector giving the Kriging variance at each test location.}
  \item{K}{ If \code{fullCov=TRUE}, the covariance matrix between the test points.}
}
\details{
This function is only for the binary LOOCV weighting method as the variance can only be obtained for this method. Mean prediction using other weighting methods can be obtained with the function \code{pred_comb_mean}.
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
