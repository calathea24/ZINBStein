#' @import GeneNet
covcal <- function(data) {
  x = wt.scale(as.matrix(data))
  n = nrow(x)
  sw = sqrt(rep(1/n, n))
  S = crossprod(sweep(x, MARGIN=1, STATS=sw, FUN="*"))
  return(S)
}

#' @import GeneNet
GeneNet <- function(data) {
  return(estimate.lambda(data, verbose = FALSE))
}

Fisher2011 <- function(data) {
  x = as.matrix(data)
  n = nrow(x)
  p = ncol(x)

  S.emp = covcal(x)
  #Shrinkage intensity estimate (equations 11)
  a1 = sum(diag(S.emp))/p
  a2 = (n^2/((n-1)*(n+2)*p))*(sum(diag(S.emp %*% S.emp)) - 1/n*(sum(diag(S.emp))^2))

  lambda = (1/n*a2 + p/n*(a1^2))/((n+1)/n*a2+p/n*(a1^2)-2*a1+1)
  lambda = max(0, min(1, lambda)) ## truncate lambda to [0,1]
  return(lambda)
}

HimenoYamada2014 = function(data) {
  x = as.matrix(data)
  N = nrow(x)
  n = N - 1
  p = ncol(x)

  #Shrinkage intensity estimate based on theorem 1 from T.Himeno, T.Yamada 2014 and function of Lambda shrinking
  #towards identity matrix (page 254, A.Touloumis, 2015)
  S = covcal(x)
  Q = 1/(N-1)*sum((diag(x%*%t(x)))^2)
  const = (N-1)/(N*(N-2)*(N-3))
  trSigma2 = const*((N-1)*(N-2)*sum(diag(S%*%S)) + (sum(diag(S)))^2 - N*Q)
  tr2Sigma = const*(2*sum(diag(S%*%S)) +(n^2 - 3*n + 1)*(sum(diag(S)))^2 - n*Q)
  Y2N = trSigma2
  Y1N = sum(diag(S))
  Y1N.2 = tr2Sigma
  beta2 = Y2N + Y1N.2
  delta2 = N*Y2N + Y1N.2 - (N-1)*(2*Y1N - p)
  lambda = beta2/delta2
  lambda = max(0, min(1, lambda))
  return(lambda)
}

ShrinkIKS <- function(data) {
  x = as.matrix(data)
  N = nrow(x)
  n = N - 1
  p = ncol(x)

  #Shrinkage intensity estimate (equations 11)
  S = covcal(x)
  Q = 1/(N-1)*sum((diag(x%*%t(x)))^2)
  const = (N-1)/(N*(N-2)*(N-3))

  a2c = const*((N-1)*(N-2)*sum(diag(S%*%S)) + (sum(diag(S)))^2 - N*Q)/p
  a1.2 = (sum(diag(S))/p)^2

  numerator = (sum(diag(S%*%S)))/p - a2c
  denominator = (sum(diag(S%*%S)))/p - a1.2

  lambda = numerator/denominator
  lambda = max(0, min(1, lambda))
  return(lambda)
}

#' @import GeneNet
SteinShrink <- function(data, method, lambda = NULL, lambda.return = FALSE) {
  if (is.null(lambda)) {
    shrinkageFunction <- match.fun(method)
    lambda <- shrinkageFunction(data)
  }

  if (lambda.return) {
    return(lambda)
  } else {
    x = as.matrix(data)
    p = ncol(x)
    S = covcal(x)
    TS = diag(p)
    S.est = (1-lambda)*S + lambda*TS
    return(cor2pcor(S.est))
  }
}

#' @import huge
LassoShrink <- function(data, lambda = seq(1,0,length=200), method.lasso = "glasso"){
  return(huge(data, lambda = lambda, method = method.lasso, scr = FALSE, verbose = FALSE))
}



