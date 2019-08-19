#' Convert BED to BigWig
#'
#' @importFrom data.table fread
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges coverage
#' @importFrom GenomicRanges tile
#' @importFrom IRanges IRanges
#' @importFrom IRanges Views
#' @importFrom IRanges viewSums
#' @importFrom S4Vectors mcols
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom rtracklayer export.bw
#'
#' @param bed BED file path, which will be converted to BigWig
#' @param chromSizes GRanges store chromosome size
#' @param binSize bin size
#' @param minCell fiter cluster with cell less than minCell
#' @param norm BPM, readInPeak for samples of varying quality
#' @param unionPeaks unin peak sets in GRanges
#' @param scale scale factor
#' @param ... other params
#' @rdname bedToBigWig
#' @export
bedToBigWig <- function(bed, chromSizes = NULL,
                        binSize = 100,
                        minCell = 30,
                        norm = "BPM",
                        unionPeaks = NULL,
                        scale = 1e7,
                        verbose = FALSE, ...){

  if (is.null(chromSizes)){
    stop("chromSize is required")
  }

  if ( norm == "readInPeak" && is.null(unionPeaks)){
    stop("unionPeaks is required")

  }

  df <- fread(bed, header = FALSE, data.table = FALSE)

  gr <- GRanges(
    seqnames = df[,1],
    IRanges(df[,2]+1, df[,3])
  )

  cvg <- coverage(gr)
  windows <- unlist(tile(chromSizes, width = binSize))
  cov_by_wnd <- Views(cvg, windows)
  sum_by_wnd <- viewSums(cov_by_wnd)

  if (norm == "readInPeak"){
    # calculate size factor
    sizeFactor <- readInPeakNorm(gr, unionPeaks = unionPeaks,
                                 scale = scale)
    if (verbose){
      message(sprintf("%s of %s", bed, sizeFactor))
    }
    mcols(windows)['score'] <- unlist(sum_by_wnd) * sizeFactor
  } else{
    mcols(windows)['score'] <- unlist(sum_by_wnd) / sum(sum(cvg)) * scale
  }

  seqlengths(windows) <- end(chromSizes)

  bwOut <- gsub("bed","bw",bed)
  export.bw(windows, bwOut)
}


#' Calculate the factor of read in peak normalization
#'
#' @importFrom GenomicRanges tile
#' @importFrom S4Vectors mcols
#' @param fragmentFile file path of Tn5 offset-corrected(+4 for +strand and -5 for -strand) insertion sites
#'
#' @param unionPeaks union peak
#' @param chromSize chrome size
#'
#' @rdname readInpeakNorm
#' @export
readInPeakNorm <- function(gr, unionPeaks,
                           scale = 1e7,...){


  # Count the insertion in peak ---------------------------------------------
  mcols(gr)['RG'] <- "pseduoBulk"
  matPeaks <- countInsertions(unionPeaks, gr, 'RG')[[1]]


  # calculate size factor ---------------------------------------------------
  readsInpeaks <- sum(Matrix::colSums(matPeaks))
  # Normalize the total number of reads by a scale factor that
  # converts all samples to a constant N million reads within peaks
  scaleFactor <- scale / readsInpeaks

  return(scaleFactor)

}

