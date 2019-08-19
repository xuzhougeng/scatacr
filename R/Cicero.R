#' GRange to Feature
#'
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start
#' @importFrom BiocGenerics end
#'
#' @rdname grToFeature
#' @export
grToFeature <- function(gr){
  peakinfo <- data.frame(
    row.names = paste(seqnames(gr),start(gr),end(gr),sep="_"),
    site_name = paste(seqnames(gr),start(gr),end(gr),sep="_"),
    chr = gsub("chr","",as.character(seqnames(gr))),
    bp1 = start(gr),
    bp2 = end(gr)
  )
  return(peakinfo)
}


#' feature to GRanges
#'
#' @importFrom stringr str_split
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#'
#' @rdname featureToGR
#' @export
featureToGR <- function(feature){
  featureSplit <- str_split(paste0(feature),
                            pattern = "_", n = 3, simplify = TRUE)
  gr <- GRanges(featureSplit[,1],
                IRanges(as.integer(featureSplit[,2]),
                        as.integer(featureSplit[,3])))
  return(gr)
}

#' make cell_data_set object from SummarizedExperiment object
#'
#' @importFrom monocle3 new_cell_data_set
#' @importFrom monocle3 fData
#' @param se SummarizedExperiment object.
#' @param binarize binary the matrix.
#'
#' @rdname makeCDS
#' @export
makeCDS <- function(se, binarize = TRUE, ...){

  peakinfo <- grToFeature(se)
  mat <- assay(se)
  if(binarize){
    mat@x[which(mat@x > 0)] <- 1
  }
  cellinfo <- data.frame(colData(se))
  cellinfo$cells <- rownames(cellinfo)
  cds <-  suppressWarnings(new_cell_data_set(
    mat,
    cellinfo,
    peakinfo
  ))
  fData(cds)$chr <- as.character(fData(cds)$chr)
  fData(cds)$bp1 <- as.numeric(as.character(fData(cds)$bp1))
  fData(cds)$bp2 <- as.numeric(as.character(fData(cds)$bp2))
  cds <- cds[order(fData(cds)$chr, fData(cds)$bp1),]
  return(cds)

}

#' make Cicero cell_data_set object from cell_data_setobject
#'
#' @importFrom monocle3 detect_genes
#' @importFrom monocle3 estimate_size_factors
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment reducedDims
#' @importFrom cicero make_cicero_cds
#' @param se SummarizedExperiment object.
#' @param k Number of cells to aggregate per bin.
#'
#' @rdname makeCicero
#' @export
makeCicero <- function(cds, k = 50, ...){
  # detect gene
  cds <- detect_genes(cds)
  cds <- estimate_size_factors(cds)

  #Reduced Dimensions
  dimred <- data.frame(
    row.names = colnames(cds),
    colData(cds)$UMAP1,
    colData(cds)$UMAP2
  )
  reducedDims(cds)$UMAP <- dimred[colnames(cds),]

  ciceroObj <- make_cicero_cds(cds,
                               k = k,
                               reduced_coordinates = dimred[colnames(cds),])


  return(ciceroObj)
}


#' Computing grouped correlations
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment rowData
#' @importFrom GenomicRanges resize
#' @importFrom GenomicRanges findOverlaps
#' @param ciceroObj cicero cell_data_set object
#' @rdname calcConnections
#' @export
calcConnections <- function(ciceroObj,
                            flank = 250*10^3,
                            ...
                            ){

  gr <- featureToGR(rowData(ciceroObj)[[1]])

  o <- suppressWarnings(as.matrix(
    findOverlaps(resize(resize(gr,1,"center"), 2*flank + 1, "center"),
                 resize(gr,1,"center"),
                 ignore.strand=TRUE) ))

  o <- data.table::as.data.table(
    data.frame(i = matrixStats::rowMins(o),
               j = matrixStats::rowMaxs(o)))
  o <- data.frame(o[!duplicated(o),])
  o <- o[o[,1]!=o[,2],]

  o$cor <- rowCorCpp(o[,1], o[,2], assay(ciceroObj), assay(ciceroObj))

  connections <- data.frame(
    Peak1 = rowData(ciceroObj)[[1]][o[,1]],
    Peak2 = rowData(ciceroObj)[[1]][o[,2]],
    coaccess = o[,3]
  )
  return(connections)
}


#' get TxDb genes
#'
#' @importFrom IRanges subsetByOverlaps
#' @importFrom AnnotationDbi mapIds
#' @importFrom S4Vectors mcols
#' @param txdb txdb object
#' @param orgdb orgdb object
#' @param gr GRanges
#' @param ignore.strand ignore the strand
#'
#' @rdname getTxDbenes
#' @export
getTxDbGenes <- function(txdb = NULL,
                         orgdb = NULL,
                         gr = NULL,
                         ignore.strand = TRUE){

  if (is.null(genome)) {
    if (is.null(txdb) | is.null(orgdb)) {
      stop("If no provided genome then you need txdb and orgdb!")
    }
  }

  if (is.null(gr)) {
    genes <- GenomicFeatures::genes(txdb)
  }else {
    genes <- suppressWarnings(
      subsetByOverlaps(GenomicFeatures::genes(txdb),
                       gr,
                       ignore.strand = ignore.strand))
  }

  if (length(genes) > 1) {
    mcols(genes)$symbol <- suppressMessages(mapIds(orgdb,
                                                   keys = mcols(genes)$gene_id,
                                                   column = "SYMBOL",
                                                   keytype = "ENTREZID",
                                                   multiVals = "first"))
    genes <- sort(sortSeqlevels(genes), ignore.strand = TRUE)
    names(genes) <- NULL
    out <- genes
  }else {
    out <- GRanges(seqnames(gr),
                   ranges = IRanges(0, 0),
                   gene_id = 0,
                   symbol = "none")[-1]
  }

  return(out)

}


