#' @import scran
scaling.fac <- function(data, clusters = TRUE, bulk = FALSE){
  #Details: p x n count data, estimate factor for each cell/observation
  #         bulk method - divide library sizes to mean
  #         scalingfNoClustes - no cluster specified
  #         scalingfClusters - calculate within and between clusters
  if (bulk){
    lib.sizes <- colSums(data)
    scales <- lib.sizes/mean(lib.sizes)
  } else {
    sce <- SingleCellExperiment(assays = list(counts = data))
    if (clusters){
      clusters <- suppressWarnings(quickCluster(sce, method = "hclust", min.size = 10))
      sce <- computeSumFactors(sce, clusters = clusters)
    } else {
      sce <- computeSumFactors(sce)
    }
    scales <- sizeFactors(sce)
  }
  return(scales)
}

copulaNB <- function(xs, scales) {
  xs <- as.numeric(xs)
  nblike <- function(ps) {
    -sum(dnbinom(xs, mu=ps[1]*scales, size=ps[2], log = T))
  }
  ml <- optim(c(mean(xs),1), nblike, method = "L-BFGS-B", lower = c(0.001,0.001), upper = c(Inf,Inf))
  mu <- ml$par[1]
  si <- ml$par[2]

  p <- pnbinom(xs, size = si, mu = mu)
  p[p==0] = 0.00001
  p[p==1] = 0.99999
  co <- qnorm(p)
  return(co)
}

copulaECDF <- function(xs) {
  emp.cdf <- ecdf(xs)
  p <- emp.cdf(xs)
  p[p==0] = 0.00001
  p[p==1] = 0.99999
  co <- qnorm(p)
  return(co)
}

#' @import huge
preprocessing <- function(data, method, clusters = FALSE, bulk = TRUE) {
  x <- t(data)

  if (method %in% c("scalingfCluster", "scalingfNoCluster", "scalingfBulk")) {
    clusters <- method %in% c("scalingfCluster", "scalingfBulk")
    bulk <- method %in% c("scalingfBulk")
    weights <- scaling.fac(x, clusters = clusters, bulk = bulk)
    xs <- t(sweep(x, 2, weights, "/"))
    return(xs)

  } else if (method == "logT") {
    return(log(data + 1))

  } else if (method %in% c("copulaNB", "copulaECDF")) {
    NB <- method == "copulaNB"
    weights <- scaling.fac(x, clusters = FALSE, bulk = FALSE)
    da <- t(do.call(rbind, lapply(1:nrow(x), function(i) {
      if (NB) {
        copulaNB(x[i,], weights)
      } else {
        copulaECDF(x[i,])
      }
    })))
    rownames(da) <- rownames(x)
    colnames(da) <- colnames(x)
    return(da)

  } else if (method == "nonparanormal") {
    data.npn <- huge.npn(data, npn.func = "truncation", verbose = FALSE)
    return(data.npn)

  } else if (method == "none") {
    return(data)

  }
}














