#' @import GeneNet
# Significance test of pcor --------------------
ptDistribution <- function(pcor) {
  pvalue <- fdrtool(pcor, statistic = "correlation", verbose = FALSE, plot = FALSE)$pval
  return(pvalue)
}

#' @import GeneNet
#' @export
pshrunkv1 <- function(pcor, p, n, lambda) {
  # https://doi.org/10.1093/bioinformatics/btz357
  #=====================Distribution function=====================#
  F0 <- function(x, k) {
    fp = ((k-3)/2)*log((1-lambda)^2-x^2) - log(beta(1/2,(k-1)/2)) - (k-2)*log(1-lambda)
    return(fp)
  }

  #=====================Estimate k with MLE=====================#
  ### Simulate data with all partial correlation = 0
  ### Fitting distribution function on it

  k.est <- function(p, n, lambda) {
    # Simulate data
    pcor.sim = ggm.simulate.pcor(p, 0)
    data.sim = ggm.simulate.data(n, pcor.sim)
    pcor.est = pcor.shrink(data.sim, lambda, verbose = FALSE)
    pcor.shrunk = sm2vec(pcor.est)

    #MLE
    F0.mle <- function(k0) {
      -sum(F0(pcor.shrunk, k = k0))
    }
    ml <- optim(100, F0.mle, method="L-BFGS-B", lower=5, upper=Inf)

    return(ml$par[1])
  }

  ### Iteration 50 to take mean value of k in each simulation
  k.iter <- 50
  k <- mean(sapply(1:k.iter, function(i){k.est(p, n, lambda)}))

  #=====================Calculate p-value=====================#
  pval.mat = matrix(Inf, length(pcor), 1)
  fp0 = function(x) {(((1-lambda)^2-x^2)^((k-3)/2))/(beta(1/2,(k-1)/2)*((1-lambda)^(k-2)))}
  for (i in 1:length(pcor)) {
    int = integrate(fp0, lower = -(1-lambda), upper = -abs(pcor[i]))
    pval.mat[i] = 2*int$value
  }
  return(pval.mat)
}

#' @import GeneNet
#' @import GeneNetTools
#' @export
pshrunkv2 <- function(pcor, p, n, lambda) { # Codes from GeneNetTools package

  k = k.shrunk(p, n, lambda) ## estimate by standard

  # Estimate p values using the shrunk t-test

  # Rescale pcor
  rr <-  unlist(pcor) / (1- lambda)

  # T-test
  tt <- rr * sqrt( (k-1)/(1 - rr^2 ) )

  # p-vals student
  pval.shrunk <- 2*pt(q = -abs(tt),
                      df = (k - 1),
                      lower.tail = TRUE, log.p = FALSE)
  pval.shrunk <- as.matrix(pval.shrunk)
  return(pval.shrunk)
}

#' @import GeneNet
#' @export
pmontecarlo <- function(pcor, p, n, graph, model, preProcess, lambda, ShrinkMet, algorithm, number = 50, depth = 50000) { ## generalize: required data simulation model and pcor.est model

  cum.pv<-matrix(0,length(pcor),1)

  # Simulate null hypothetic GGM coefficients for "number" times
  for (i in 1:number){
    sim <- data.simulate(p, n, graph, degree = 0, type = model, depth = depth)
    r.data <- sim$sim.dat
    pcor.data <- sim$pcor.sim

    if (model == "ESCO"){
      methods = strsplit(preProcess,"-")[[1]]
      for (i in 1:length(methods)){
        met = methods[[i]]
        if (met == "scalingfCluster" | met == "scalingfNoCluster" | met == "scalingfBulk"){
          r.data = preprocessing(r.data, method = met)
        } else {
          r.data = r.data[,colnames(pcor.data)]
          r.data = preprocessing(r.data, method = met)
        }
      }
      r.data = r.data[,colnames(pcor.data)]
    }

    if (ShrinkMet == "Stein") {
      r.monte.GGM = SteinShrink(r.data, method = algorithm, lambda = lambda)
    } else {
      message("pMonteCarlo only for Stein-type estimators")
    }
    r.monte <- sm2vec(r.monte.GGM)

    # compare the real coefficients against r.monte
    pv <- sapply(pcor, function(x) sum(abs(r.monte)>=abs(x))/length(r.monte))
    cum.pv <- cum.pv + pv
  }

  # p values
  p.monte<-cum.pv/number
  return (p.monte)
}

#' @import GeneNet
#' @export
pcor.testing <- function(pc, ..., pmodel = "ptDistribution", method = "BH")
{
  cutoff.pval <- 0.05
  p <- ncol(pc)
  pcor <- sm2vec(pc)

  ##1/Computing p-values
  pCalculation <- match.fun(pmodel)
  argumentList <- list(pcor = pcor, ...)
  argumentRequire <- names(formals(pCalculation))
  pval <- do.call(pCalculation, argumentList[argumentRequire])

  ##2/FDR method
  if (method == "localFDR"){
    p.adj <- fdrtool(pcor, statistic = "correlation", verbose = FALSE, plot = FALSE)$lfdr  #posterior probabilities (= 1- local fdr)
  } else if (method == "tailFDR"){
    p.adj <- fdrtool(pcor, statistic = "correlation", verbose = FALSE, plot = FALSE)$qval
  } else {
    p.adj <- p.adjust(pval, method = method)
  }

  ##3/ Returning adjacency matrix
  indexes = sm.index(pc)
  colnames(indexes) = c("node1", "node2")
  result = cbind(indexes, p.adj)
  idx = which(p.adj <= cutoff.pval)
  adj.mat = matrix(0, nrow = p, ncol = p)
  for (i in idx) {
    node1 = result[i,1]
    node2 = result[i,2]
    adj.mat[node1,node2] = 1
  }
  adj.mat = adj.mat + t(adj.mat)
  return(adj.mat)
}



