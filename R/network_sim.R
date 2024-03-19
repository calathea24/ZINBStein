#' Simulate scale-free network following Barabasi-Albert model
#'
#' Simulate scale-free network using Barabasi-Albert model from https://doi.org/10.1103/RevModPhys.74.47
#' with two important properties. Growth and preferential attachment
#'
#' @param num_of_nodes integers. Number of nodes or features in network
#' @param m can be floating numbers, lambda value of Poisson distribution which simulates k, number of iterations for preferential attachments
#'
#' @return adjacency matrix
#'
#' @export
# Output of network simulation: adjacency matrix, no self-connection

simulateBAmodel <- function(num_of_nodes, m) {
  N <- num_of_nodes

  # 1. Generate base graph with 2 nodes and 1 edge between them
  adj_mat <- matrix(0, nrow = N, ncol = N)
  adj_mat[2,1] <- adj_mat[1,2] <- 1
  degs <- c(1,2)

  # 2.  Add nodes gradually until total nodes equal N
  for (i in 3:N) {
    k <- rpois(n = 1, lambda = m)
    alist <- vector()
    for (j in 1:k) {
      a <- sample(degs, size = 1)
      adj_mat[i,a] <- adj_mat[a,i] <- 1
      alist <- append(alist, c(a,i))
    }
    degs <- append(degs, alist) ## preferential attachment
  }
  return(adj_mat)
}



#' Network simulation
#'
#' Simulating different network structure (random, hub, cluster, band, scale-free). Simulation of hub, cluster and band follow models from huge.generator (huge R package)
#'
#' @param d number of nodes in network
#' @param graph type of network structure. Five options: random (default), hub, cluster, band, scale-free
#' @param degree average number of edges per node in random graph. Default: 2
#' @param g number of hubs or clusters in hub and cluster graphs, respectively.
#' @param prob parameter for cluster graph, probability that a pair of nodes has an edge in each cluster. Default: 6*g/d
#' @param pcor.return TRUE or FALSE. return only partial correlation matrix or a list of all relevant simulated matrices.
#'
#' @return return partial correlation matrix by default
#'
#' @example
#'
#' @export
# Network structure simulation ------------------
graph.generator <- function(d = NULL, graph = "random", degree = 2, g = NULL, prob = NULL, pcor.return = TRUE)
{
  # modify to sync three parameters, degree, etaA and prob
  gcinfo(FALSE)

  if(is.null(g)){ ## group in hub and cluster graph
    g = 1
    if(graph == "hub" || graph == "cluster"){
      if(d > 40)  g = ceiling(d/20)
      if(d <= 40) g = 2
    }
  }

  if(graph == "cluster"){
    if(is.null(prob)){
      if(d/g > 30)  prob = 0.3
      if(d/g <= 30)  prob = min(1,6*g/d)
    }
    prob = sqrt(prob/2)*(prob<0.5)+(1-sqrt(0.5-0.5*prob))*(prob>=0.5)
  }

  ####################
  eps = 0.0001
  ###################

  # parition variables into groups
  g.large = d%%g
  g.small = g - g.large
  n.small = floor(d/g)
  n.large = n.small+1
  g.list = c(rep(n.small,g.small),rep(n.large,g.large))
  g.ind = rep(c(1:g),g.list)
  rm(g.large,g.small,n.small,n.large,g.list)
  gc()


  # Build the graph structure with adjacency matrix
  adj.mat = matrix(0,d,d);

  if(graph == "random") {
    num.edges = d*(d-1)/2
    degree = degree
    etaA = (degree*d)/num.edges
    num.elements = ceiling(num.edges*etaA)
    element.idx = sample(1:num.edges, num.elements)
    adj.lo = rep(0,num.edges)
    adj.lo[element.idx] = 1
    adj.mat[lower.tri(adj.mat)] = adj.lo
    adj.mat = adj.mat + t(adj.mat)
  }

  if(graph == "band"){
    for(i in 1:g){
      diag(adj.mat[1:(d-i),(1+i):d]) = 1
      diag(adj.mat[(1+i):d,1:(d-1)]) = 1
    }
  }

  if(graph == "cluster"){
    for(i in 1:g){
      tmp = which(g.ind==i)
      tmp2 = matrix(runif(length(tmp)^2,0,0.5),length(tmp),length(tmp))
      tmp2 = tmp2 + t(tmp2)
      adj.mat[tmp,tmp][tmp2<prob] = 1
      rm(tmp,tmp2)
      gc()
    }
  }

  if(graph == "hub"){
    for(i in 1:g){
      tmp = which(g.ind==i)
      adj.mat[tmp[1],tmp] = 1
      adj.mat[tmp,tmp[1]] = 1
      rm(tmp)
      gc()
    }
  }

  if(graph == "scale-free"){
    adj.mat = simulateBAmodel(d, 1.5)
  }

  diag(adj.mat) = 0

  # Simulate precision matrix from adjacency matrix
  selected.edges = which(adj.mat[upper.tri(adj.mat)]!=0)
  precision = matrix(0, nrow=d, ncol=d)
  precision[upper.tri(precision)][selected.edges] = runif(length(selected.edges),-1.0,+1.0)
  precision = precision + t(precision)

  #Constructing diagonally dominant/ positive definite matrix
  for(i in 1:d)
  {
    diag(precision)[i] = sum(abs(precision[,i])) + eps
  }

  #Converting to partial correlation matrix
  pcor = cov2cor(precision) # Standardize precision matrix
  # change signs of the off-diagonal entries to obtain pcor matrix
  pcor = -pcor
  diag(pcor) = -diag(pcor) # keep positive sign on diagonal

  if(pcor.return){
    return(pcor)
  } else {
    m = pcor
    m = -m
    diag(m) = -diag(m)
    cov = pseudoinverse(m)
    return(list(adj.matrix = adj.mat, precision.matrix = precision, pcor.matrix = pcor, cov.matrix = cov))
  }
}








