#' Run MACS2 to call peak in pseduo-bulk
#'
#' @param bed.path input BED path
#' @param peak.dir output peak dir
#' @param macs2Path macs2 path
#' @param method q qvalue, p pvalue
#' @param cutoff cutoff of pvalue or qvalue
#' @param shift The arbitrary shift in bp
#' @param extsize The arbitrary extension size in bp
#' @param genomeSize Effective genome size
#' @param bedGraph output bedGraph
#' @param addtion additional options passed to MACS2
#'
#' @export
runMACS2 <- function(bed.path = NULL,
                     peak.dir = NULL,
                     macs2Path = "macs2",
                     method = "q",
                     cutoff = 0.05,
                     shift  =  -75,
                     extsize = 150,
                     bedGraph = TRUE,
                     genomeSize = 2.7e9,
                     addtion = NULL){

  if ( is.null(bed.path)){
    stop("bed.path can not be empty")
  }

  if ( is.null(peak.dir)){
    stop("peak.dir can not be empty")
  }

  if ( ! dir.exists(peak.dir)){
    dir.create(peak.dir, recursive = TRUE)
  }

  cmdPeaks <- sprintf(
    "%s callpeak -g %s --name %s --treatment %s --outdir %s --format BED --nomodel --call-summits --nolambda --keep-dup all",
    macs2Path,
    genomeSize,
    sub("\\.bed","", basename(bed.path)),
    bed.path,
    peak.dir
  )
  if (!is.null(shift) & !is.null(extsize)) {
    cmdPeaks <- sprintf("%s --shift %s --extsize %s", cmdPeaks, shift, extsize)
  }
  if (tolower(method) == "p") {
    cmdPeaks <- sprintf("%s -p %s", cmdPeaks, cutoff)
  }else {
    cmdPeaks <- sprintf("%s -q %s", cmdPeaks, cutoff)
  }

  if (bedGraph){
    cmdPeaks <- sprintf("%s -B --SPMR", cmdPeaks)
  }

  if ( ! is.null(addtion) ){
    cmdPeaks <- sprintf("%s %s", cmdPeaks, addtion)
  }

  message("Running Macs2...")
  message(cmdPeaks)
  system(cmdPeaks, intern = TRUE)

}


# remove the peak which overlap with the more significant peak
#' only keep the most significant peak in overlap region
#'
#' @importFrom GenomicRanges reduce
#' @importFrom GenomicRanges findOverlaps
#' @importFrom GenomicRanges reduce
#' @importFrom GenomeInfoDb sortSeqlevels
#' @importFrom IRanges subsetByOverlaps
#'
#' @param gr GRanges object
#' @param by the overlapped peak order standard
#' @param decreasing logical. Should the sort order be increasing or decreasing?
#' @param verbose  Be chatty and report processing?
#'
#' @export
nonOverlappingGRanges <- function(gr, by = "score",
                                  decreasing = TRUE,
                                  verbose = FALSE){
  stopifnot(by %in% colnames(mcols(gr)))

  # cluster the peak
  clusterGRanges <- function(gr,
                             filter = TRUE,
                             by = "score",
                             decreasing = TRUE, ...){
    gr <- sort(sortSeqlevels(gr))
    r <- reduce(gr, min.gapwidth=0L, ignore.strand=TRUE)
    o <- findOverlaps(gr,r)
    mcols(gr)$cluster <- subjectHits(o)
    gr <- gr[order(mcols(gr)[,by], decreasing = decreasing),]
    gr <- gr[!duplicated(mcols(gr)$cluster),]
    gr <- sort(sortSeqlevels(gr))
    mcols(gr)$cluster <- NULL
    return(gr)
  }


  if(verbose){
    message("Converging", appendLF = FALSE)
  }
  i <-  0
  gr_converge <- gr
  while(length(gr_converge) > 0){
    if(verbose){
      message(".", appendLF = FALSE)
    }
    i <-  i + 1
    gr_selected <- clusterGRanges(gr = gr_converge,
                                  filter = TRUE,
                                  by = by,
                                  decreasing = decreasing)
    #blacklist selected gr
    gr_converge <- subsetByOverlaps(gr_converge ,
                                    gr_selected, invert=TRUE)

    #if i=1 then set gr_all to clustered
    if(i == 1){
      gr_all <- gr_selected
    }else{
      gr_all <- c(gr_all, gr_selected)
    }
  }
  if(verbose){
    message("\nSelected ", length(gr_all), " from ", length(gr))
  }
  gr_all <- sort(sortSeqlevels(gr_all))
  return(gr_all)
}



#' extent summit to build a peak set
#'
#' @importFrom IRanges subsetByOverlaps
#' @importFrom GenomicRanges GRangesList
#' @importFrom GenomicRanges resize
#'
#' @param summits summit file generated from MACS2 callpeak
#' @param nSummits number of summit to use
#' @param extend extend the summit
#' @param chromSizes chromosome length
#' @param blacklist genome black list in BED format
#'
#' @export
summit2peak <- function(summits,
                        nSummits = 200000,
                        extend = 250,
                        chromSizes = NULL,
                        blacklist = NULL){

  # check file
  stopifnot(all(sapply(summits, file.exists)))

  # Deal with blacklist -----------------------------------------------------
  if(is.null(blacklist)){
    blacklist <- GRanges()
  }else if(is.character(blacklist)){
    blacklist <- import.bed(blacklist)
  }

  gr_list <- lapply(summits, function(x){
    summit <- readSummits(x)
    # extent the summit
    summit <- resize(summit, width = 2 * extend + 1, fix = "center" )
    # remove the peak beyond the end of chromosome
    summit <- subsetByOverlaps(summit, chromSizes, type="within")
    # remove the peak overlap with blacklist region
    summit <- subsetByOverlaps(summit, blacklist,invert=TRUE)
    # only keep the most significant peak in the clustered peak set
    submit <- nonOverlappingGRanges(summit, by="score", decreasing=TRUE)

    # order the summit by the score
    submit <- submit[order(submit$score,decreasing=TRUE)]

    if(!is.null(nSummits)){
      submit <- head(submit, nSummits)
    }

    # peak score normalization:
    # quantile = trunc(rank(v))/length(v)
    # v represents the vector of MACS2 peaks scores
    peak_scores <- mcols(submit)$score
    quantile <- trunc(rank(peak_scores))/length(peak_scores)
    mcols(submit)$scoreQuantile <- quantile

    return(submit)
  })

  gr_list <- GRangesList(gr_list)

  #Non Overlapping of each cluster
  grNonOverlapping <- nonOverlappingGRanges(unlist(gr_list),
                                            by = "scoreQuantile",
                                            decreasing = TRUE)

  grNonOverlapping
}


#' merge the peak from all sample
#'
#' @importFrom GenomeInfoDb sortSeqlevels
#' @importFrom GenomicRanges GRangesList
#'
#' @param ... peak from different sample
#' @param by the overlapped peak order standard
#' @param decreasing logical. Should the sort order be increasing or decreasing?
#' @param verbose  Be chatty and report processing?
#'
#' @export
mergePeak <- function(...,
                      by = "scoreQuantile",
                      decreasing = TRUE){

  gr_list <- GRangesList(...)

  gr_final <- nonOverlappingGRanges(unlist(gr_list),
                                   by = by,
                                   decreasing = decreasing)

  gr_final <- sort(sortSeqlevels(gr_final))
  return(gr_final)

}


