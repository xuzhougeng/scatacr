#' Write the fragments per cluster
#'
#' @importFrom BiocGenerics start
#' @importFrom BiocGenerics end
#' @importFrom GenomeInfoDb seqnames
#' @importFrom SummarizedExperiment colData
#' @importFrom magrittr %>%
#'
#' @param x Seurat or SummarizedExperiment object
#' @param fragmentFiles vector store path of fragments GRanges Rds
#' @param dirClusters output directory
#'
#' @rdname writeFragmentsToBed
#' @export
writeFragmentsToBed <- function(x,
                                fragmentFiles = NULL,
                                dirClusters = "results/LSI-Cluster-Beds/"){

  # check -------------------------------------------------------------------
  if (is.null(fragmentFiles)){
    stop("fragmentFiles can not be empty")
  }

  if ( ! dir.exists(dirClusters)){
    dir.create(dirClusters, recursive = TRUE)
  }

  # get the cell of each cluster  -------------------------------------------
  if ( inherits(x, "Seurat") ){
    clusterResults <- split(rownames(x@meta.data),
                            paste0("Cluster",
                                   x@meta.data[,ncol(x@meta.data)]))
  } else if( inherits(x, "SummarizedExperiment" ) ) {
    clusterResults <- split(rownames(colData(x)),
                            paste0("Cluster",colData(x)[,"Clusters"]))
  } else{
    stop("x is not Seurat or SummarizedExperiment object")
  }

  for(i in seq_along(fragmentFiles)){
    fragments <- readRDS(fragmentFiles[i])
    for(j in seq_along(clusterResults)){
      message(sprintf("%s of %s", j, length(clusterResults)))
      fragmentsj <- fragments[fragments$RG %in% clusterResults[[j]]]
      if(length(fragmentsj) > 0){
        out <- data.frame(
          chr = c(seqnames(fragmentsj), seqnames(fragmentsj)),
          # BED start with 0
          start = c(as.integer(start(fragmentsj) - 1), as.integer(end(fragmentsj) - 1)),
          end = c(as.integer(start(fragmentsj)), as.integer(end(fragmentsj)))
        ) %>% readr::write_tsv(
          x = .,
          append = TRUE,
          path = file.path(dirClusters, paste0(names(clusterResults)[j], ".bed")),
          col_names = FALSE)
      }
    }
  }
}
