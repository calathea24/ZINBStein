#' Simulate homogeneous scRNAseq data from ESCO model
#'
#' Simulate a simple model of single-cell RNA sequencing (scRNAseq) data
#' with cells from the same population and no outlier cells (homogeneous) based on ESCO model (ESCO R package)
#'
#' @param GeneNumbers integer, number of genes to be simulated
#' @param CellNumbers integer, number of homogeneous cells to be simulated
#' @param pcorSim partial correlation matrix which contains information of networks. Non-zero value equals to presence of edge between two genes.
#' @param depth integer, sequencing depth. For scRNAseq, 20 000 reads are recommended by most workflows. Default: 50 000 reads
#' @param zeroprop value from 0 to 1, maximum proportion of zero counts of genes in network (pcorSim). Default: 0.4
#' @param zisim TRUE or FALSE. Simulating zero-inflated counts from model of ESCO R package
#' @param verbose TRUE or FALSE. Information of simulating process
#'
#' @return list of simulated count (genes x cells) and partial correlation matrix.
#'
#' @examples
#'
#'
#' @import ESCO
#' @import foreach
#' @import GeneNet

escoSimulateSingleGroup <- function(GeneNumbers, CellNumbers, pcorSim, depth = 50000, zeroprop = 0.4, zisim = FALSE, verbose = FALSE)
{
  params <- newescoParams()
  checkmate::assertClass(params, "escoParams")
  params <- setParams(params,
                      nGenes = GeneNumbers,
                      nCells = CellNumbers,
                      lib.loc = round(log(as.numeric(depth)), digits = 1),
                      lib.scale = 0.1,
                      out.prob = 0)
  validObject(params)
  nCells <- getParam(params, "nCells")
  nGenes <- getParam(params, "nGenes")

  #1-Set up Simulation Object
  # Set up name vectors
  cell_names <- paste0("Cell", seq_len(nCells))
  gene_names <- paste0("Gene", seq_len(nGenes))
  # Create SingleCellExperiment to store simulation
  cells <-  data.frame(Cell = cell_names)
  rownames(cells) <- cell_names
  features <- data.frame(Gene = gene_names)
  rownames(features) <- gene_names
  sim <- SingleCellExperiment(rowData = features, colData = cells, metadata = list(Params = params))

  #2-Simulate library sizes
  sim <- ESCO::escoSimLib(sim, verbose)

  #3-Simulate base gene means
  sim <- ESCO::escoSimGeneMeans(sim, verbose)

  #4-Simulate gene means
  sim <- ESCO::escoSimSingleCellMeans(sim, verbose)

  #5-Simulate true counts
  # bcv (Biological Coefficient of Variation) simulation
  bcv_common <- getParam(params, "bcv.common")
  bcv_df <- getParam(params, "bcv.df")
  basecell_means <- assays(sim,withDimnames = FALSE)$BaseCellMeans
  basecell_means <- as.matrix(basecell_means)

  bcv <- (bcv_common + (1 / sqrt(basecell_means)))*sqrt(bcv_df/rchisq(nGenes, df = bcv_df))
  dimnames(bcv) <- dimnames(basecell_means)
  bcv <- as.matrix(bcv)

  cell_means <- matrix(rgamma(nGenes * nCells, shape = 1 / (bcv ^ 2),
                              scale = basecell_means * (bcv ^ 2)),
                       nrow = nGenes, ncol = nCells)

  # correlated gene simulation, change randcor function in ESCO to simulate from partial correlation matrix
  randcor <- function(pcorr){
    corrgenes <- sample(1:nGenes, nrow(pcorSim))
    corr <- pcor2cor(pcorr)
    colnames(corr) = rownames(corr) = gene_names[corrgenes]
    colnames(pcorr) = rownames(pcorr) = gene_names[corrgenes]
    re = list("pcor" = pcorr, "cor" = corr, "corrgenes" = corrgenes)
    return(re)
  }

  corr_sim <- randcor(pcorSim)
  rho <- corr_sim[["cor"]]
  pcor.truth <- corr_sim[["pcor"]]
  corrgenes <- corr_sim[["corrgenes"]]
  params <- setParams(params, corr = list(rho))
  copular <- ESCO:::randcop(rho, nCells)

  mulpo <- function(i) {
    count = rnbinom(nGenes, size = 1/(bcv[,i]^2), mu = basecell_means[,i])
    count[corrgenes] = qnbinom(copular[,i],  size = 1/(bcv[corrgenes,i]^2), mu = basecell_means[corrgenes,i])
    return(count)
  }

  # change ESCO code to do sequential processing instead of parallel, reduce error when simulating
  total <- ncol(basecell_means)
  if (verbose) {
    pb <- progress_bar$new(
      format = "progress = :letter [:bar] :elapsed | eta: :eta", total = total, width = 60)
    progress <- function(n){
      pb$tick(tokens = list(letter = rep("", total)[n]))
    }
    opts <- list(progress = progress)
    true_count <- foreach(i = 1:total, .combine = cbind, .options.snow = opts, .export = c("rowData")) %do% {
      return(mulpo(i))
    }
  } else{
    true_count <- foreach(i = 1:total, .combine = cbind, .export = c("rowData")) %do% {
      return(mulpo(i))
    }
  }
  colnames(true_count) <- cell_names
  rownames(true_count) <- gene_names
  true_count[true_count==Inf] <- max(true_count[true_count!=Inf]) + 10

  # make sure correlated genes having zero counts less than certain percent
  if(!is.null(zeroprop)){
    data_check <- t(true_count[corrgenes,])
    zero_percent <- sapply(1:ncol(data_check), function(x){
      length(which(as.numeric(data_check[,x])==0))/nrow(data_check)
    })
    names(zero_percent) <- colnames(data_check)
    idx <- colnames(data_check)[which(zero_percent >= zeroprop)]
    for (i in idx){
      z_value <- zero_percent[i]
      increment <- 0
      while (z_value >= zeroprop) {
        increment <- increment + 1
        size <- 1/(bcv[i,]^2) + increment
        mu <- basecell_means[i,] + increment
        true_count[i,] <- qnbinom(copular[i,], size = size , mu = mu )
        z_value <- length(which(as.numeric(true_count[i,]) == 0))/ncol(true_count)
      }
    }
  }

  assays(sim,withDimnames = FALSE)$TrueCounts <- true_count
  assays(sim,withDimnames = FALSE)$CellMeans <- cell_means
  metadata(sim)$Params = params

  #6-Simulate zero inflation
  if (zisim) {
    dropout_mid <- rep(getParam(params, "dropout.mid"), nCells)
    dropout_shape <- rep(getParam(params, "dropout.shape"), nCells)
    dropout_cort <- getParam(params, "dropout.cort")
    cell_normmeans <- median(colSums(cell_means))*t(t(cell_means)/colSums(cell_means))
    dropout_shape <- getParam(params, "dropout.shape")

    # Generate probabilites based on expression
    logistic <- function(x, x0, k) {
      1 / (1 + exp(-k * (x - x0)))
    }
    drop_prob <- vapply(seq_len(nCells), function(idx) {
      eta <- log(cell_normmeans[,idx])
      return(logistic(eta, x0 = dropout_mid[idx], k = dropout_shape[idx]))
    }, FUN.VALUE = rep(0,nGenes))

    if(!dropout_cort) keep_prob <- 1 - drop_prob
    keep_prob[keep_prob>1] <- 1
    keep_prob[keep_prob<0] <- 0
    keep_prob[is.na(keep_prob)] <- 1
    keep <- matrix(rbinom(nCells * nGenes, 1, keep_prob),
                   nrow = nGenes, ncol = nCells)
    ZI_counts <- true_count * keep

    return(list(TrueCounts = true_count, ZICounts = ZI_counts, pcor.truth = pcor.truth))
  } else {
    return(list(TrueCounts = true_count, pcor.truth = pcor.truth))
  }
}



#' Simulate normally distributed or scRNAseq data based on partial correlation matrix
#'
#' Using partial correlation matrix to simulate data
#'
#' @param features integer. Number of genes or features
#' @param samples integer. Number of cells/samples/observations
#' @param graph type of graph or network. Available options: random, hub, cluster, band, scale-free
#' @param degree degree of simulated graph or average number of edges per node. Default: 2
#' @param type normal or ESCO for normally distributed or scRNAseq data simulation, respectively
#' @param ZI TRUE or FALSE (default). Whether simulating zero-inflated count or not in scRNAseq data simulation
#' @param depth integer, sequencing depth. For scRNAseq, 20 000 reads are recommended by most workflows. Default: 50 000 reads
#'
#' @return list of simulated partial correlation matrix (pcor.sim) and data matrix (sim.dat)
#'
#' @examples
#'
#' @import GeneNet
#' @export

data.simulate <- function(features, samples, graph = "random", degree = 2, type = "normal", ZI = FALSE, depth = 50000)
{
  ## Step 1: simulate partial correlation matrix
  pcor.sim = graph.generator(d = features, graph = graph, degree = degree, pcor.return = TRUE)

  ## Step 2: simulate data using pcor.sim
  if (type == "normal") {
    sim.dat = ggm.simulate.data(sample.size = samples, pcor = pcor.sim)
  } else if (type == "ESCO") {
    sd.check = TRUE
    while (sd.check){ ## Make sure genes in network are not all zeros or having 0 standard deviation
      sim = escoSimulateSingleGroup(15000, samples, pcor.sim, zisim = ZI, depth = depth)
      sim.dat = t(sim$TrueCounts)
      if (ZI) {
        sim.dat = t(sim$ZICounts)
      }
      pcor.sim = sim$pcor.truth
      rm(sim)
      gc()
      dat = sim.dat[,colnames(pcor.sim)]
      sd.vec = apply(dat,2,sd)
      sd.check = any(sd.vec==0)
    }
  } else {
    message("Current simulation does not support other models.")
  }
  return(list(pcor.sim = pcor.sim, sim.dat = sim.dat))
}










