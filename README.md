# combinedkriging : A R package for Bayesian optimization using a combination of Kriging models with random length-scales

## Summary :

`combinedkriging` is a R package that was developed to perform Bayesian Optimization (BO) of high-dimensional black-box functions.
Classical Efficient Global Optimization (EGO) using Kriging models typically struggles in high dimensions because of the difficulties in estimating the length-scales hyperparameters (using maximum likelihood estimaiton for example).

`combinedkriging` provides a method to build a combination of Kriging sub-models, with random length-scales, which can be used for high-dimensional BO. In particular, it gives functions to randomly sample appropriate length-scales for the sub-models, to build the combination with various weighting methods, and to perform BO using the combination. It also includes trust region implementations to improve the performances of BO.

## Installation :

`combinedkriging` can be installed from GitHub, for the very latest version:

``` r
# If not already installed, install package `devtools` with `install.packages("devtools")`
devtools::install_github("TAppriou/combinedkriging")
```
It requires the following dependencies:

- `DiceKriging` : available in the CRAN
- `copula` : available in the CRAN
``` r
install.package(c("DiceKriging","copula"))
```
- `nestedKriging` : available at https://github.com/drulliere/nestedKriging

## Simple use example

``` r

library(combinedkriging)

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

```


