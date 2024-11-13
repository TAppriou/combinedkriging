\name{run_EGO}
\alias{run_EGO}
\title{EGO algorithm}
\description{
Run the EGO algorithm
}
\usage{
run_EGO(n_steps, fun, ..., X_train, Y_train, model, bornes_theta_opt=c(0.1,20), estim_mean=TRUE, combination=TRUE, n_tree_levels=NULL, optimizer="cmaes", TR_type="TREGO", nugget_OK=1e-6, nugget_comb=0, log_file=NULL)
}
\arguments{
  \item{n_steps}{Number of EGO iterations.}
  \item{fun}{A function to be minimized, with first argument the vector of parameters over which minimization is to take place. It should return a scalar result.}
  \item{\dots}{Further arguments to be passed to \code{fun}.}
  \item{X_train}{The training points. The ith row contains the d input variables corresponding to the ith observation.}
  \item{Y_train}{The function values at the training points. The ith element of the vector is the function value for the ith training point.}
  \item{model}{A combination of Kriging sub-models obtained with the function \code{combined_krg}, or an ordinary Kriging model obtained with \code{krg_MLE}.}
  \item{bornes_theta_opt}{Bounds used for the length-scale MLE optimization.}
  \item{estim_mean}{Do we use a non-zero constant trend? Default value is \code{TRUE}}.
  \item{combination}{\code{TRUE} for a combined Kriging model, \code{FALSE} for an ordinary Kriging model.}
  \item{n_tree_levels}{Number of levels in the binary tree for the combination. The corresponding number of sub-models is \code{2^(n_tree_levels - 1)}.}
  \item{optimizer}{The optimizer for the EI. Only "cmaes" is implemented in this version.}
  \item{TR_type}{The type of Trust region used. Options are no trust region "\code{NoTR}", TREGO implementation of trust regions "\code{TREGO}", TURBO implementation of trust regions "\code{TURBO}".}
  \item{nugget_OK}{Optional nugget added to the diagonal of the correlation matrix.}
  \item{log_file}{An optional folder to write EGO progress. If \code{log_file=NULL}, the progress is written in the console instead.}
}
\value{
A list containing:
  \item{model}{The combined Kriging or ordinary Kriging model at the end of the EGO iterations,}
  \item{X}{All design points after the EGO iterations. The ith row contains the d input variables corresponding to the ith point.}
  \item{Y}{All function values at the design points after the EGO iterations. The ith element of the vector is the function value for the ith point.}
  \item{time}{A list giving the overall time for building the models and for optimizing the EI.}
  \item{history_center_TR}{A matrix giving the history for the trust region center. The ith row contains the d parameters for the TR center at the ith iteration,}
  \item{history_length_TR}{A vector giving the length of the trust region at each EGO iteration.}
}

\examples{

### An example for the 10D sphere function ###

# Test function: Sphere function
fun_sphere <- function(x) {
  return(sqrt(sum((x-0.5)^2)))
}

d <- 15 #Dimension
n_train_init <- 20 #Number of training samples
n_steps <- 30
n_tree_levels <- 5 #Number of levels for the binary tree combinations
n_comb_max <- 2^(n_tree_levels - 1)
nugget <- 1e-4

set.seed(123) #Set the random seed

# Random training points
X_train_init <- matrix(0,nrow=n_train_init,ncol=d)
for (i in 1:n_train_init) {
  X_train_init[i,] <- runif(d, min=0, max=1)
}

# Training values
Y_train_init <- apply(X_train_init, MARGIN=1, fun_sphere)


# EGO for the combination

time1 <- Sys.time()

# Lower and upper bounds for sampling the length-scales
bounds_comb <- bornes_theta(d=d, X=X_train_init, kernel="matern5_2")

# Sampling the length-scales
theta_list <- sample_theta_entropy(X=X_train_init, n_comb=n_comb_max, bornes=bounds_comb,
                                   kernel="matern5_2", iso=TRUE)

# Build the initial combination
model_comb <- combined_krg(X_train=X_train_init, Y_train=Y_train_init, theta_list=theta_list,
                           kernel="matern5_2", estim_mean=TRUE, nugget=0,
                           n_tree_levels=n_tree_levels, weighting_method="binLOO")

# Run EGO
EGO_comb <- run_EGO(n_steps, fun=fun_sphere, X_train=X_train_init, Y_train=Y_train_init, model=model_comb, estim_mean=TRUE,
                    combination=TRUE, n_tree_levels=n_tree_levels, optimizer="cmaes",
                    TR_type="TREGO", nugget=0, log_file=NULL)

time2 <- Sys.time()
time_comb <- difftime(time2,time1,units="secs")

# Evolution of best value
best_sol_comb <- rep(NA, 1+n_steps)
best_sol_comb[1] <- min(Y_train_init)
for (i in 1:n_steps) {
  best_sol_comb[i+1] <- min(best_sol_comb[i], EGO_comb$Y[n_train_init+i])
}


# EGO for ordinary Kriging


time3 <- Sys.time()


# Build the initial Kriging model using MLE
model_krg <- krg_MLE(X_train=X_train_init, Y_train=Y_train_init, bornes_opt=c(0.1,20),
                     kernel="matern5_2", theta_init=NULL, estim_mean=TRUE, nugget=NULL)

# Run EGO
EGO_krg <- run_EGO(n_steps, fun=fun_sphere, X_train=X_train_init, Y_train=Y_train_init, model=model_krg, bornes_theta_opt=c(0.1,20),
                   estim_mean=TRUE, combination=FALSE, optimizer="cmaes",
                   TR_type="TREGO", nugget=NULL, log_file=NULL)

time4 <- Sys.time()
time_krg <- difftime(time4,time3,units="secs")

# Evolution of best value
best_sol_krg <- rep(NA, 1+n_steps)
best_sol_krg[1] <- min(Y_train_init)
for (i in 1:n_steps) {
  best_sol_krg[i+1] <- min(best_sol_krg[i], EGO_krg$Y[n_train_init+i])
}

}
