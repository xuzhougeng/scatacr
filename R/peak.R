#
#' Extend the peak
#'
#' @importFrom rtracklayer import.bed
#' @importFrom GenomicRanges GRanges
#' @importFrom magrittr %>%
#' @importFrom GenomicRanges resize
#' @importFrom S4Vectors mcols
#' @importFrom IRanges subsetByOverlaps
#'
#' @param df data.frame object include peak information
#' @param genomeSize genome size GRange object
#' @param extend extend the summit
#' @param blacklist genome black list
#' @param nSummits number of summit
#' @param ... Arguments passed to seurat RunLSI
#'
#' @return Seurat object
#'
#' @rdname extendedPeakSet
#' @export
extendedPeakSet <- function(df, genomeSize = NULL,
                            extend = 250, blacklist = NULL,
                            nSummits = 100000, ...){

  #Check-------
  stopifnot(extend > 0)
  stopifnot("samples" %in% colnames(df))
  stopifnot("groups" %in% colnames(df))
  stopifnot("summits" %in% colnames(df))
  stopifnot(!is.null(genomeSize))
  stopifnot(all(apply(df,1,function(x){file.exists(paste0(x[3]))})))
  #------------
  #Deal with blacklist
  if(is.null(blacklist)){
    blacklist <- GRanges()
  }else if(is.character(blacklist)){
    blacklist <- import.bed(blacklist)
  }
  stopifnot(inherits(blacklist,"GenomicRanges"))

  #------------
  #Time to do stuff
  chromSizes <- GRanges(names(genomeSize),
                        IRanges(1, genomeSize))
  chromSizes <- GenomeInfoDb::keepStandardChromosomes(chromSizes,
                                                      pruning.mode = "coarse")
  groups <- unique(df$groups)

  groupGRList <- GenomicRanges::GenomicRangesList(lapply(seq_along(groups), function(i){
    df_group = df[which(df$groups==groups[i]),]
    grList <- GenomicRanges::GenomicRangesList(lapply(paste0(df_group$summits), function(x){
      extended_summits <- readSummits(x) %>%
        resize(., width = 2 * extend + 1, fix = "center") %>%
        subsetByOverlaps(.,chromSizes,type="within") %>%
        subsetByOverlaps(.,blacklist,invert=TRUE) %>%
        nonOverlappingGRanges(., by="score", decreasing=TRUE)
      extended_summits <- extended_summits[order(extended_summits$score,decreasing=TRUE)]
      if(!is.null(nSummits)){
        extended_summits <- head(extended_summits, nSummits)
      }
      mcols(extended_summits)$scoreQuantile <- trunc(rank(mcols(extended_summits)$score))/length(mcols(extended_summits)$score)
      extended_summits
    }))

    #Non Overlapping
    grNonOverlapping <- nonOverlappingGRanges(unlist(grList),
                                              by = "scoreQuantile",
                                              decreasing = TRUE)

    #Free Up Memory
    remove(grList)
    grNonOverlapping

  }))
  grFinal <- nonOverlappingGRanges(unlist(groupGRList),
                                   by = "scoreQuantile",
                                   decreasing = TRUE)
  grFinal <- sort(sortSeqlevels(grFinal))
  return(grFinal)
}


#Helper Functions
readSummits <- function(file){
  df <- suppressMessages(data.table::fread(file,
                                           col.names = c("chr","start","end","name","score"),
                                           data.table = FALSE))
  df <- df[,c(1,2,3,5)] #do not keep name column it can make the size really large
  return(GenomicRanges::makeGRangesFromDataFrame(df=df,keep.extra.columns = TRUE,starts.in.df.are.0based = TRUE))
}


nonOverlappingGRanges <- function(gr, by = "score",
                                  decreasing = TRUE,
                                  verbose = FALSE){
  stopifnot(by %in% colnames(mcols(gr)))

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
    gr_selected <- clusterGRanges(gr = gr_converge, filter = TRUE, by = by, decreasing = decreasing)
    gr_converge <- subsetByOverlaps(gr_converge ,gr_selected, invert=TRUE) #blacklist selected gr
    if(i == 1){ #if i=1 then set gr_all to clustered
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

#' cluster GRanges
#'
#' @importFrom GenomicRanges reduce
#' @importFrom GenomicRanges findOverlaps
#' @importFrom GenomeInfoDb sortSeqlevels
#'
clusterGRanges <- function(gr,
                           filter = TRUE,
                           by = "score",
                           decreasing = TRUE, ...){
  gr <- sort(sortSeqlevels(gr))
  r <- GenomicRanges::reduce(gr, min.gapwidth=0L, ignore.strand=TRUE)
  o <- findOverlaps(gr,r)
  mcols(gr)$cluster <- subjectHits(o)
  gr <- gr[order(mcols(gr)[,by], decreasing = decreasing),]
  gr <- gr[!duplicated(mcols(gr)$cluster),]
  gr <- sort(sortSeqlevels(gr))
  mcols(gr)$cluster <- NULL
  return(gr)
}

#' Run MACS2
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
#' @rdname runMACS2
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


#' createPeakSummarizedExperiment
#'
#' @importFrom GenomeInfoDb keepSeqlevels
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics Reduce
#' @importFrom BiocGenerics start
#' @importFrom BiocGenerics end
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors SimpleList
#'
#' @importFrom magrittr %>%
#' @param fragmentFiles vector include
#' @param genome genomeSize GRanges, used in extendedPeakSet
#' @param extend extend, used in extendedPeakSet
#' @param blacklist blacklist, use in extendedPeakSet
#' @param nSummits  nSummits, use in extendedPeakSet
#' @param chromKeep chromosome to keep
#' @param clusters clusters generated by Seurat SNN graph clustering
#'
#' @rdname createPeakSummarizedExperiment
#' @export
createPeakSummarizedExperiment <- function(fragmentFiles,
                                           genome,
                                           groups = "scATAC",
                                           extend = 250,
                                           blacklist = NULL,
                                           nSummits = 200000,
                                           chromKeep = NULL,
                                           dirPeaks = "results/LSI-Cluster-Peaks",
                                           clusters = NULL, ...){

  df <- data.frame(
    samples = gsub("\\_summits.bed","",
                   list.files(dirPeaks, pattern = "\\_summits.bed", full.names = FALSE)),
    groups = groups,
    summits = list.files(dirPeaks, pattern = "\\_summits.bed", full.names = TRUE)
  )

  unionPeaks <- extendedPeakSet(
    df = df,
    genomeSize = genome,
    extend = extend,
    blacklist = blacklist,
    nSummits = nSummits
  )

  # filter chromosome like ChrY, ChrM
  if ( ! is.null(chromKeep)){
    unionPeaks <- unionPeaks[seqnames(unionPeaks) %in% chromKeep ]
    unionPeaks <- keepSeqlevels(unionPeaks,  chromKeep )
  }

  #Create Counts list
  countsPeaksList <- lapply(seq_along(fragmentFiles), function(i){
    message(sprintf("%s of %s", i, length(fragmentFiles)))
    gc()
    countInsertions(unionPeaks, readRDS(fragmentFiles[i]), by = "RG")
  })

  #CountsMatrix
  mat <- lapply(countsPeaksList, function(x) x[[1]]) %>% Reduce("cbind",.)
  frip <- lapply(countsPeaksList, function(x) x[[2]]) %>% unlist
  total <- lapply(countsPeaksList, function(x) x[[3]]) %>% unlist
  se <- SummarizedExperiment(
    assays = SimpleList(counts = mat),
    rowRanges = unionPeaks
  )
  rownames(se) <- paste(seqnames(se),start(se),end(se),sep="_")
  colData(se)$FRIP <- frip
  colData(se)$uniqueFrags <- total / 2

  if ( ! is.null(clusters)){
    colData(se)$Clusters <- clusters
  }

  return(se)

}








