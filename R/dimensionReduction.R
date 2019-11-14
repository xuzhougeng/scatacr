#' LSI reduction
#'
#' @importFrom irlba irlba
#'
#' @param mat dgCMatrix object
#' @param nComponents Number of singular values to compute
#' @param binarize binarize the matrix
#' @param nFeatures Getting top N features for  SVD
#' if nFeatures < 1, Getting the features  expressing >= cell * nFeatures
#' @param ... Arguments passed to seurat RunLSI
#'
#' @return SVD matrix
#'
#' @export
lsi <- function(mat, nComponents = 50,
                      binarize = TRUE,
                      nFeatures = 2000,
                      ...){

  #TF IDF LSI adapted from flyATAC
  cs <- Matrix::colSums(mat)
  if(binarize){
    message(paste0("Binarizing matrix..."))
    mat@x[mat@x > 0] <- 1
  }
  if(!is.null(nFeatures)){

    if (nFeatures < 1){
      numCellsCounted <-  Matrix::rowSums(mat)
      threshold <-  ncol(mat) * nFeatures
      message(paste0("Getting features open in >= ", threshold, "..."))
      mat  <-  mat[numCellsCounted >= threshold,]
    } else{
      message(paste0("Getting top ", nFeatures, " features..."))
      mat <- mat[head(order(Matrix::rowSums(mat),decreasing = TRUE),
                      nFeatures),]
    }

  }

  row.names(mat) <- seq(1, nrow(mat))

  #Calc TF IDF
  message("Computing Term Frequency IDF...")
  freqs <- Matrix::t(Matrix::t(mat)/Matrix::colSums(mat))
  idf   <- as(log(1 + ncol(mat) / Matrix::rowSums(mat)), "sparseVector")
  tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% freqs

  #Calc SVD then LSI
  message("Computing SVD using irlba...")
  svd <- irlba(tfidf,
               nComponents, nComponents,
               maxit = 1000)
  svdDiag <- matrix(0, nrow=nComponents, ncol=nComponents)
  diag(svdDiag) <- svd$d
  matSVD <- t(svdDiag %*% t(svd$v))
  rownames(matSVD) <- colnames(mat)
  colnames(matSVD) <- paste0("LSI_",seq_len(ncol(matSVD)))

  return(matSVD)


}

