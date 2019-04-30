[![Travis build status](https://travis-ci.org/arnaldpuy/sensobol.svg?branch=master)](https://travis-ci.org/arnaldpuy/sensobol) [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/arnaldpuy/sensobol?branch=master&svg=true)](https://ci.appveyor.com/project/arnaldpuy/sensobol) [![CRAN_Status_Badge_version_last_release](https://www.r-pkg.org/badges/version-last-release/sensobol)](https://cran.r-project.org/package=sensobol) [![metacran downloads](https://cranlogs.r-pkg.org/badges/last-week/sensobol)](https://cran.r-project.org/package=sensobol) [![metacran downloads](https://cranlogs.r-pkg.org/badges/grand-total/sensobol)](https://cran.r-project.org/package=sensobol) [![Coverage status](https://codecov.io/gh/arnaldpuy/sensobol/branch/master/graph/badge.svg)](https://codecov.io/github/arnaldpuy/sensobol?branch=master) [![license](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2579855.svg)](https://doi.org/10.5281/zenodo.2579855)

# sensobol

The goal of `sensobol` is to provide a set of functions to swiftly compute and visualize up to third-order Sobol' sensitivity indices. The functions allow to: 
- Create the sample matrices for the model evaluation.
- Compute and bootstrap up to third-order effects.
- Assess the approximation error of Sobol' indices.
- Plot the model uncertainty and the Sobol' indices.

## Installation
To install the stable version on [CRAN](https://CRAN.R-project.org/package=sensobol), use

```r
install.packages("sensobol")
```
To install the development version, use devtools:

``` r
install.packages("devtools") # if you have not installed devtools package already
devtools::install_github("arnaldpuy/sensobol", build_vignettes = TRUE)
```

## Example

This brief example shows how to compute Sobol' indices. For a more detailed explanation of the package functions, check the vignette.

``` r
## Load the package:
library(sensobol)

## Create sample matrix to compute first, total and second-order indices:
A <- sobol_matrices(n = 1000, k = 3,  second = TRUE)

## Compute the model output (using the Ishigami test function):
Y <- ishigami_Mapply(A)

## Compute the Sobol' indices (first, total and second-order):
sens <- sobol_indices(Y = Y, params = colnames(data.frame(A)), R = 100, n = 1000, second = TRUE)
```

## Citation

Please use the following citation if you use `sensobol` in your publications:

```r
Arnald Puy (2019). sensobol: Computation of High-Order Sobol' Sensitivity Indices. R package
  version 0.2.0 http://github.com/arnaldpuy/sensobol
```

A BibTex entry for LaTex users is:

```r
@Manual{,
    title = {sensobol: Computation of High-Order Sobol' Sensitivity Indices},
    author = {Arnald Puy},
    year = {2019},
    note = {R package version 0.2.0},
    url = {http://github.com/arnaldpuy/sensobol},
  }
```
