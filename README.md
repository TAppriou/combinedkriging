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
``` r
install.package("DiceKriging")
```
- `nestedKriging` : available at https://github.com/drulliere/nestedKriging

