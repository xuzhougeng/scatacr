#' LSI reduction analysis
#'
#' @importFrom Seurat CreateSeuratObject
#' @importFrom Seurat RunLSI
#'
#' @param mat GenomicRanges object  to store feature information
#' @param nComponents Number of singular values to compute
#' @param binarize binarize the matrix
#' @param nFeatures Getting top N features for download analysis
#' @param ... Arguments passed to seurat RunLSI
#'
#' @return Seurat object
#'
#' @rdname seuratLSI
#' @export
seuratLSI <- function(mat, nComponents = 50,
                        binarize = TRUE,
                        nFeatures = NULL, ...){

  #TF IDF LSI adapted from flyATAC
  cs <- Matrix::colSums(mat)
  if(binarize){
    message(paste0("Binarizing matrix..."))
    mat@x[mat@x > 0] <- 1
  }
  if(!is.null(nFeatures)){
    message(paste0("Getting top ", nFeatures, " features..."))
    mat <- mat[head(order(Matrix::rowSums(mat),decreasing = TRUE),nFeatures),]
  }

  # create a seurat object
  row.names(mat) <- seq(1, nrow(mat))
  obj <- CreateSeuratObject(mat,
                            assay = "ATAC",
                            project='scATAC',
                            min.cells=0, min.features=0)
  #Calc TF IDF, Calc SVD then LSI
  obj <- RunLSI(object = obj, n = nComponents, scale.max = NULL)

  return(obj)

}

#' Find Cluster using SNN graph clustering methods
#'
#' @importFrom Seurat FindNeighbors
#' @importFrom Seurat FindClusters
#'
#' @param obj Matrix  cell-by-feature matrix
#' @param minGroupSize the mininum cell required to be cluster,
#' 200 is required according to 2019, NBT.
#' Set it to NULL if the peaks are called
#' @param dims.use Dimensions of reduction to use as input
#' @param initialResolution Value of the resolution parameter,
#' use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities
#' @param n.start Number of random starts
#' @param ... Arguments passed to Seurat's FindNeighbors or FindClusters
#'
#' @return Seurat object
#'
#' @rdname addClusters
#' @export
addClusters <- function(obj,
                        minGroupSize = 50,
                        dimsUse = seq_len(50),
                        initialResolution = 0.8,
                        nStart = 10, ...){
  # get the SNN
  obj<- FindNeighbors(obj, reduction = "lsi", dims = dimsUse)

  #First Iteration of Find Clusters

  currentResolution <- initialResolution
  obj <- FindClusters(obj, resolution = currentResolution, n.start = nStart)
  res_name <-  colnames(obj@meta.data)[grepl(currentResolution, colnames(obj@meta.data))]
  minSize <- min(table(obj@meta.data[[res_name]]))
  nClust <- length(unique(paste0(obj@meta.data[[res_name]])))
  message(sprintf("Current Resolution = %s, No of Clusters = %s, Minimum Cluster Size = %s", currentResolution, nClust, minSize))

  if (is.null(minGroupSize)){
    return(obj)
  }

  #If clusters are smaller than minimum group size
  while(minSize <= minGroupSize){

    res_name <-  colnames(obj@meta.data)[grepl(currentResolution, colnames(obj@meta.data))]

    obj@meta.data <- obj@meta.data[,-which(colnames(obj@meta.data)==res_name)]
    currentResolution <- currentResolution*initialResolution
    obj <- FindClusters(obj, resolution = currentResolution)

    res_name <-  colnames(obj@meta.data)[grepl(currentResolution, colnames(obj@meta.data))]
    minSize <- min(table(obj@meta.data[[res_name]]))
    nClust <- length(unique(paste0(obj@meta.data[[res_name]])))

    message(sprintf("Current Resolution = %s, No of Clusters = %s, Minimum Cluster Size = %s", currentResolution, nClust, minSize))
  }
  return(obj)
}




# seuratLSI <- function(mat, nComponents = 50, binarize = TRUE, nFeatures = NULL){
#   #TF IDF LSI adapted from flyATAC
#   cs <- Matrix::colSums(mat)
#   if(binarize){
#     message(paste0("Binarizing matrix..."))
#     mat@x[mat@x > 0] <- 1
#   }
#   if(!is.null(nFeatures)){
#     message(paste0("Getting top ", nFeatures, " features..."))
#     mat <- mat[head(order(Matrix::rowSums(mat),decreasing = TRUE),nFeatures),]
#   }
#
#   #Calc TF IDF
#   message("Computing Term Frequency IDF...")
#   freqs <- t(t(mat)/Matrix::colSums(mat))
#   idf   <- as(log(1 + ncol(mat) / Matrix::rowSums(mat)), "sparseVector")
#   tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% freqs
#   #Calc SVD then LSI
#   message("Computing SVD using irlba...")
#   svd <- irlba::irlba(tfidf, nComponents, nComponents)
#   svdDiag <- matrix(0, nrow=nComponents, ncol=nComponents)
#   diag(svdDiag) <- svd$d
#   matSVD <- t(svdDiag %*% t(svd$v))
#   rownames(matSVD) <- colnames(mat)
#   colnames(matSVD) <- paste0("PC",seq_len(ncol(matSVD)))
#   #Make Seurat Object
#   message("Making Seurat Object...")
#   mat <- mat[1:100,] + 1
#   rownames(mat) <- seq(1,nrow(mat))
#   obj <- CreateSeuratObject(mat, project='scATAC', min.cells=0, min.features=0)
#   obj <- CreateDimReducObject(object = obj, reduction.type = "pca", slot = "cell.embeddings", new.data = matSVD)
#   obj <- CreateDimReducObject(object = obj, reduction.type = "pca", slot = "key", new.data = "PC")
#
#   obj@reductions$pca <- CreateDimReducObject(embeddings = matSVD, key = "PC", assay = "scATAC")
#   return(obj)
# }

# addClusters <- function(obj, minGroupSize = 50, dims.use = seq_len(50), initialResolution = 0.8){
#   #First Iteration of Find Clusters
#   currentResolution <- initialResolution
#   obj <- FindClusters(object = obj, reduction.type = "pca", dims.use = dims.use, resolution = currentResolution, print.output = FALSE)
#   minSize <- min(table(obj@meta.data[[paste0("res.",currentResolution)]]))
#   nClust <- length(unique(paste0(obj@meta.data[[paste0("res.",currentResolution)]])))
#   message(sprintf("Current Resolution = %s, No of Clusters = %s, Minimum Cluster Size = %s", currentResolution, nClust, minSize))
#   #If clusters are smaller than minimum group size
#   while(minSize <= minGroupSize){
#     obj@meta.data <- obj@meta.data[,-which(colnames(obj@meta.data)==paste0("res.",currentResolution))]
#     currentResolution <- currentResolution*initialResolution
#     obj <- FindClusters(object = obj, reduction.type = "pca", dims.use = dims.use, resolution = currentResolution, print.output = FALSE, force.recalc = TRUE)
#     minSize <- min(table(obj@meta.data[[paste0("res.",currentResolution)]]))
#     nClust <- length(unique(paste0(obj@meta.data[[paste0("res.",currentResolution)]])))
#     message(sprintf("Current Resolution = %s, No of Clusters = %s, Minimum Cluster Size = %s", currentResolution, nClust, minSize))
#   }
#   return(obj)
# }
