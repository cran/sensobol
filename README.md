[![Travis build status](https://travis-ci.org/arnaldpuy/sensobol.svg?branch=master)](https://travis-ci.org/arnaldpuy/sensobol) [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/arnaldpuy/sensobol?branch=master&svg=true)](https://ci.appveyor.com/project/arnaldpuy/sensobol) [![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/sensobol)](https://cran.r-project.org/package=sensobol) [![license](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2579855.svg)](https://doi.org/10.5281/zenodo.2579855)

# sensobol

The goal of `sensobol` is to provide a set of functions to swiftly compute and visualize up to third-order Sobol' sensitivity indices. The functions allow to: 
- Create the sample matrices for the model evaluation.
- Compute and bootstrap up to third-order effects.
- Assess the approximation error of Sobol' indices.
- Plot the model uncertainty and the Sobol' indices.

## Installation

You can install the released version of sensobol as follows:

``` r
install.packages("devtools") # if you have not installed devtools package already
devtools::install_github("arnaldpuy/sensobol", build_vignettes = TRUE)
```

## Example

This brief example shows how to compute Sobol' indices. For a more detailed explanation of the package functions, check the [vignette](https://github.com/arnaldpuy/sensobol/blob/master/vignettes/sensobol.Rmd).

``` r
## Create sample matrix
A <- sobol_matrices(n = 1000, 
k = 3, 
second = TRUE)

## Compute the model output (using the Ishigami test function):
Y <- ishigami_Mapply(A)

## Compute the Sobol' indices:
sens <- sobol_indices(Y = Y, 
params = colnames(data.frame(A)),
R = 100, 
n = 1000, 
second = TRUE)
```

