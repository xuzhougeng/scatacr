#' Build window along the genome
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges tile
#' @importFrom GenomeInfoDb keepStandardChromosomes
#'
#' @param seqname sequence name
#' @param size sequence size
#' @param windowSize each window size
#'
#' @export
buildWindows <- function(seqname,
                         size,
                         windowSize = 5000){

  chromSizes <- GRanges(seqname, IRanges(1,size))
  chromSizes <- keepStandardChromosomes(chromSizes,
                                        pruning.mode = "coarse")

  windows <- unlist(tile(chromSizes, width = windowSize))
  return(windows)
}


#' Read summit
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom data.table fread
#'
#' @param file summit file
#'
#' @export
readSummits <- function(file){
  df <- fread(file,
              col.names = c("chr","start","end","name","score"),
              data.table = FALSE,
              verbose = FALSE)
  #do not keep name column it can make the size really large
  df <- df[,c(1,2,3,5)]
  gr <- makeGRangesFromDataFrame(df=df,
                                 keep.extra.columns = TRUE,
                                 starts.in.df.are.0based = TRUE)
  return(gr)
}



#' Convert dgCMatrix into  matrix
#'
#' @param mat dgcMatrix object
#' @export
as_matrix <- function(mat){

  row_pos <- mat@i
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])

  tmp <- asMatrix(rp = row_pos, cp = col_pos, z = mat@x,
                  nrows =  mat@Dim[1], ncols = mat@Dim[2])

  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)

}

#' Write the fragments to BED
#'
#' @importFrom BiocGenerics start
#' @importFrom BiocGenerics end
#' @importFrom BiocGenerics subset
#' @importFrom GenomeInfoDb seqnames
#' @importFrom SummarizedExperiment colData
#'
#' @param gr GRanges which store the framgments postion information
#' @param subset subset the GRanges
#' @param filename output filename
#' @param append If TRUE, the file is opened in append mode and column names (header row) are not written.
#'
#' @export
writeFragmentsToBed <- function(gr,
                                subset = NULL,
                                filename =  NULL,
                                append = FALSE){

  # check -------------------------------------------------------------------
  if (is.null(filename)){
    stop("filename can not be empty")
  }

  if (!is.null(subset)) {
    gr <-  subset(gr, RG %in% subset)
  }
  out <- data.frame(
          chr = c(seqnames(gr), seqnames(gr)),
          # BED start with 0
          start = c(as.integer(start(gr) - 1), as.integer(end(gr) - 1)),
          end = c(as.integer(start(gr)), as.integer(end(gr)))
        )
  data.table::fwrite(
          x = out,
          file = filename,
          append = append,
          sep = "\t",
          col.names = FALSE)
}


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


