# Zero-inflated negative binomial modelling -----
kroneckerDelta <- function(xs){
  xs <- as.numeric(xs)
  re <- rep(0,length(xs))
  re[xs==0] <-1
  re
}

logsum <- function(a, b){
  m <- pmax(a,b)
  log(exp(a-m)+exp(b-m))+m
}

maxlikeNB <- function(xs, scaling) {
  xs <- as.numeric(xs)
  if(is.null(scaling)){
    scaling = rep(1, length(xs))
  }

  nblike <- function(ps)
  {
    w <- ps[3]
    -sum(logsum(log(1-w)+dnbinom(xs,size=ps[1],mu=ps[2]*scaling,log = T),log(w)+log(kroneckerDelta(xs))))
  }
  ml <- optim(c(1,mean(xs),0.5),nblike,method="L-BFGS-B",lower=c(0.001,0.001,0.001),upper=c(Inf,Inf,0.999))
  return(ml$par)
}

get_mix_parameters_new = function(count, scales){
  count = as.matrix(count)
  parslist = lapply(1:nrow(count), function(ii) {
    paramt = maxlikeNB(count[ii,],scaling = scales)
    return(paramt)
  })
  parslist = Reduce(rbind, parslist)
  colnames(parslist) = c("size", "mu", "rate")
  return(parslist)
}

calculate_weight_new = function (x, paramt, scales){
  pz1 = paramt[3] * kroneckerDelta(x)
  if (is.null(scales)) {
    scales = rep(1, length(x))
  }
  pz2 = (1 - paramt[3]) * dnbinom(x,size=paramt[1],mu=paramt[2]*scales,log = F)
  pz = pz1/(pz1 + pz2)
  pz[pz1 == 0] = 0
  return(cbind(pz, 1 - pz))
}

if_dropout_scimpute_new = function(mat, dthre = 0.5, scales){
  pa = get_mix_parameters_new(t(mat),scales)
  pa.check = cbind(pa, "mean" = colMeans(mat))
  pa.check = data.frame(pa.check)
  pa[,3][pa.check$rate < 0.01 & pa.check$mean < 1] = 0.9
  I = ncol(mat)
  J = nrow(mat)
  droprate = sapply(1:I, function(i) {
    if(is.na(pa[i,1])) return(rep(0,J))
    wt = calculate_weight_new(mat[,i], pa[i,],scales)
    return(wt[, 1])
  })
  dropind = 1* (droprate > dthre)
  colnames(dropind) = colnames(mat)
  return(dropind)  #List of 0 (selected genes) and 1 (excluded genes)
}




#' @import GeneNet
#' @export
count_ZI <- function(x, scaling = NULL, preProcess = NULL){
  x = as.matrix(x)
  dropout_id = if_dropout_scimpute_new(x, dthre = 0.5, scales = scaling)

  if (!is.null(scaling)){
    x = sweep(x, 1, scaling, "/")
  }

  if(!is.null(preProcess)){
    methods = strsplit(preProcess,"-")[[1]]
    for (i in 1:length(methods)){
      met = methods[[i]]
      x = preprocessing(x, method = met)
    }
  }

  x = wt.scale(x)
  x[dropout_id==1] = 0
  return(x)
}
