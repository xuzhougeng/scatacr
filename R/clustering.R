#' Find Cluster using SNN graph clustering methods
#'
#' @importFrom Seurat FindNeighbors
#' @importFrom Seurat FindClusters
#'
#' @param obj Matrix  cell-by-feature matrix
#' @param min.group.size the mininum cell required to be cluster,
#' 200 is required according to 2019, NBT.
#' Set it to NULL if the peaks are called
#' @param dims.use Dimensions of reduction to use as input
#' @param init.resolution Value of the resolution parameter,
#' use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities
#' @param n.start Number of random starts
#' @param k.param Defines k for the k-nearest neighbor algorithm
#' @param prune.SNN Sets the cutoff for acceptable Jaccard index
#' when computing the neighborhood overlap for the
#' SNN construction. Any edges with values less than
#' or equal to this will be set to 0 and removed from
#' the SNN graph. Essentially sets the strigency of pruning
#' (0 — no pruning, 1 — prune everything).
#' @param annoy.metric Distance metric for annoy. Options include:
#' euclidean, cosine, manhattan, and hamming
#' @param ... Arguments passed to Seurat's
#' FindNeighbors and FindClusters
#'
#' @return named vector
#'
#' @export
clusterBySNN <- function(obj,
                        min.group.size = 200,
                        dims.use = 1:50,
                        init.resolution = 0.8,
                        n.start = 10,
                        k.param = 20,
                        prune.SNN = 1/15,
                        annoy.metric = "euclidean",
                        ...){

  # dimenstion to use for SNN
  obj <- obj[, dims.use]
  # get the SNN
  obj <- FindNeighbors(obj, k.param=k.param,
                       prune.SNN = prune.SNN,
                       annoy.metric = annoy.metric,
                       ...)

  obj <- obj[['snn']]
  #First Iteration of Find Clusters
  currentResolution <- init.resolution

  cluster <- FindClusters(obj,
                      resolution = currentResolution,
                      n.start = n.start)

  # get
  clusterNumber <- length(levels(cluster[,1]))
  minSize <- min(table(cluster[,1]))
  message(sprintf("Current Resolution = %s, No of Clusters = %s, Minimum Cluster Size = %s",
                  currentResolution,
                  clusterNumber,
                  minSize))

  if (is.null(min.group.size)){
    cls <- cluster[,1]
    names(cls) <- row.names(cluster)
    return(cls)
  }

  #If clusters are smaller than minimum group size
  while(minSize <= min.group.size){

    currentResolution <- currentResolution * init.resolution
    cluster <- FindClusters(obj, resolution = currentResolution)

    clusterNumber <- length(levels(cluster[,1]))
    minSize <- min(table(cluster[,1]))

    message(sprintf("Current Resolution = %s, No of Clusters = %s, Minimum Cluster Size = %s",
                    currentResolution,
                    clusterNumber,
                    minSize))
  }

  cls <- cluster[,1]
  names(cls) <- row.names(cluster)
  return(cls)
  return(cls)
}
