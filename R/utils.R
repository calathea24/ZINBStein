#' @import mltools
mcc.pcor <- function(est, pcor.truth) {
  pcor.est = est[upper.tri(est)]
  pcor.truth = pcor.truth[upper.tri(pcor.truth)]
  TN = length(which(pcor.est[which(pcor.truth==0)]==0))
  TP = length(which(pcor.est[which(pcor.truth!=0)]!=0))
  FP = length(which(pcor.est[which(pcor.truth==0)]!=0))
  FN = length(which(pcor.est[which(pcor.truth!=0)]==0))
  mcc = mcc(TP=TP,TN=TN,FP=FP,FN=FN)
  return(mcc)
}

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

