
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ZINBStein

<!-- badges: start -->
<!-- badges: end -->

## Overview

Zero-inflated negative binomial modelling Stein-type shrinkage
(ZINBStein) is for estimating inverse covariance matrix from single-cell
RNA sequencing (scRNAseq) data with some properties:

-   Working in high-dimensional data (number of genes is larger than
    number of cells)
-   Ensuring final estimated covariance matrix is positive-definite
-   Using zero-inflated negative binomial modelling to stratify
    “dropout” zero counts from “true” zero counts

## Installation

You can install the development version of ZINBStein like so:

``` r
# 1st version is released and installed from gitHub
install.packages("devtools")
library(devtools)
install_github("calathea24/ZINBGraphicalModel") ## will change name of repository later on
```

## Usage

Using ZINB modelling to stratify count data

``` r
library(ZINBStein)
## basic example code
```

Using Stein-type shrinkage to estimate covariance matrix and return
adjacency matrix

``` r
library(ZINBStein)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
