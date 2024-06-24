# ZINBStein

<!-- badges: start -->

<!-- badges: end -->

ZINBStein is an R package to implement the sparse inverse covariance matrix estimation framework for single-cell RNA sequencing (scRNAseq) network analysis using zero-inflated negative bionomial modelling. Reproducible snakemake workflow is in https://github.com/calathea24/ZINBGraphicalModel

## Installation

You can install the development version of ZINBStein from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("calathea24/ZINBStein")
```

## Example


``` r
library(ZINBStein)
## Simulating data with 2000 features/genes and 100 samples/cells with ESCO simulation and 50 000 read depth
p <- 2000
n <- 100
sim <- data.simulate(features = p, samples = n, type = "ESCO", depth = 50000)
pcor.sim <- sim$pcor.sim
data.sim <- sim$sim.dat


## Normalizing scRNAseq simulated data with scaling factors
data.sim <- preprocessing(data.sim, method = "scalingfNoCluster")
# Transforming data using nonparanormal method
data.sim <- data.sim[,colnames(pcor.sim)]
data.sim <- preprocessing(data.sim, method = "nonparanormal")


## Using Stein-type shrinkage Fisher2011 algorithm to shrink sample covariance matrix and estimate partial correlation matrix
pc <- SteinShrink(x, method = "Fisher2011") # return partial correlation matrix

```
