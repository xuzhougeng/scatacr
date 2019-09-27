#' Reading Fragment Files
#'
#' @importFrom data.table fread
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb sortSeqlevels
#'
#' @param file fragments file path
#' @param minFrags filter cell with few fragments
#'
#' @rdname readFragments
#' @export
readFragments <- function(file = NULL,
                          minFrags = 100, ...){

  message("Reading in Fragment files ...")
  fragments <- fread(file, header = FALSE, data.table = FALSE)
  fragments <- GRanges(
    seqnames = fragments[,1],
    IRanges(fragments[,2]+1, fragments[,3]),
    RG = fragments[,4],
    N = fragments[,5]
  )

  if ( minFrags > 0){
    message("Filtering Lowly Represented Cells...")
    tabRG <- table(fragments$RG)
    keep <- names(tabRG)[which(tabRG >= minFrags)]
    fragments <- fragments[fragments$RG %in% keep,]
    fragments <- sort(sortSeqlevels(fragments))
    return(fragments)

  } else{
    return(fragments)
  }
}

#' Get the cell number
#' @param fragments GRanges object
#' @param minFrags filter cell with few fragments
#'
#' @rdname getCellNum
#' @export
getCellNum <- function(fragments, minFrags = 0){

  cell_number <- sum(table(fragments$RG) > minFrags)

  return(cell_number)

}

#' Filter Fragments by fragments and enrichment score
#'
#' @importFrom S4Vectors mcols
#'
#' @param fragments GRanges object
#' @param tssSingles tssSingle From getTssProfile
#' @param filterFrags threshold of fragments per cell
#' @param filterTSS threshold of enrichment score of TSS
#' @param prefix prefix add to barcode
#'
#' @rdname filterFragments
#' @export
filterFragments <- function(fragments,
                            tssSingles,
                            filterFrags = 1000,
                            filterTSS = 8,
                            prefix = NULL
                            ){

  tssSingles$cellCall <- 0
  tssSingles$cellCall[tssSingles$uniqueFrags >= filterFrags &
                        tssSingles$enrichment >= filterTSS] <- 1


  fragments <- fragments[mcols(fragments)$RG %in% rownames(tssSingles)[tssSingles$cellCall==1]]

  # Add prefix to barcode --------------------------
  if ( is.null(prefix)){
    return(fragments)
  }else{
    fragments$RG <- paste0(name,"#",fragments$RG)
    return(fragments)
  }

}



#' TSS Profile
#'
#' @importFrom magrittr %>%
#' @importFrom GenomicRanges resize
#' @importFrom GenomicFeatures transcripts
#'
#' @param fragments a GRanges object
#' @param txdb txdb object
#' @param batchSize split the cell into different batch for performance
#' @rdname getTssProfile
#' @export
getTssProfile <- function(fragments,
                          txdb,
                          batchSize = 1000){
  feature <- txdb %>%
    transcripts(.) %>%
    resize(., width = 1, fix = "start") %>%
    unique

  tssProfile <- insertionProfileSingles(feature = feature,
                                        fragments = fragments,
                                        getInsertions = TRUE,
                                        batchSize = batchSize)
  tssSingles <- tssProfile$dfall
  tssSingles$uniqueFrags <- 0
  tabRG <- table(fragments$RG)
  tssSingles[names(tabRG),"uniqueFrags"] <- tabRG
  return(tssSingles)


}



#' Generate insertion profile from fragments
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors DataFrame
#' @importFrom S4Vectors mcols
#' @importFrom Matrix sparseMatrix
#' @importFrom utils txtProgressBar
#' @importFrom magrittr %>%
#' @param fragments GenomicRanges object
#' @param feature GenomicRanges object  to store feature information
#' @param by metadata in mcols
#' @param flank flank
#' @param nrom nrom
#' @param smooth smooth
#' @param range range
#' @param batchSize Split the cell into N batch for performance
#' @param ... Arguments passed to other methods
#'
#' @return An list with sparseMatrix, FRIP and total fragments
#'
#' @rdname insertionProfileSingles
#' @export
insertionProfileSingles <- function(fragments, feature, by = "RG",
                                    getInsertions = TRUE,
                                    fix = "center",
                                    flank = 2000,
                                    norm = 100,
                                    smooth = 51,
                                    range = 100,
                                    batchSize = 100, ...){

  uniqueTags <- as.character(unique(mcols(fragments)[,by]))
  splitTags <- split(uniqueTags,
                     ceiling(seq_along(uniqueTags) / batchSize))

  pb <- txtProgressBar(min = 0, max = 100, initial = 0, style = 3)
  batchTSS <- lapply(seq_along(splitTags), function(x){
    setTxtProgressBar(pb, round(x * 100 / length(splitTags), 0))
    profilex <- insertionProfileSingles_helper(
      feature=feature,
      fragments=fragments[which(mcols(fragments)[,by] %in% splitTags[[x]])],
      by = by,
      getInsertions = getInsertions,
      fix = fix,
      flank = flank,
      norm = norm,
      smooth = smooth,
      range = range
    )

    return(profilex)
  })
  df <- lapply(batchTSS, function(x) x$df) %>% Reduce("rbind",.)
  dfall <- lapply(batchTSS, function(x) x$dfall) %>% Reduce("rbind",.)
  profileMat <- lapply(batchTSS, function(x) x$profileMat) %>% Reduce("cbind",.)
  profileMatSmooth <- lapply(batchTSS, function(x) x$profileMatSmooth) %>% Reduce("cbind",.)
  return(list(df = df, dfall = dfall, profileMat = profileMat, profileMatSmooth = profileMatSmooth))
}

#' insertionProfileSingles helper
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors DataFrame
#' @importFrom S4Vectors mcols
#' @importFrom Matrix sparseMatrix
#' @importFrom utils txtProgressBar
#' @importFrom BiocGenerics start
#' @importFrom BiocGenerics which
#'
#' @param feature GenomicRanges object  to store feature information
#' @param fragments GenomicRanges object
#' @param by metadata in mcols
#' @param flank the set flank size of TSS to find the overlap of
#' insertion
#' @param nrom nrom
#' @param smooth smooth
#' @param range range
#' @param batchSize batch size
#' @param ... Arguments passed to other methods
#'
#' @return An list with sparseMatrix, FRIP and total fragments
insertionProfileSingles_helper <- function(feature,
                                           fragments,
                                           by = "RG",
                                           getInsertions = TRUE,
                                           fix = "center",
                                           flank = 2000,
                                           norm = 100,
                                           smooth = 51,
                                           range = 100,
                                           batchSize = 100, ...){
  #Convert To Insertion Sites
  if(getInsertions){
    insertions <- c(
      GRanges(seqnames = seqnames(fragments), ranges = IRanges(start(fragments), start(fragments)), RG = mcols(fragments)[,by]),
      GRanges(seqnames = seqnames(fragments), ranges = IRanges(end(fragments), end(fragments)), RG = mcols(fragments)[,by])
    )
    by <- "RG"
  }else{
    insertions <- fragments
  }
  remove(fragments)

  #center the feature
  center <- unique(resize(feature, width = 1, fix = fix, ignore.strand = FALSE))

  #get overlaps between the feature and insertions only up to flank bp
  overlap <- DataFrame(findOverlaps(query = center,
                                    subject = insertions,
                                    maxgap = flank,
                                    ignore.strand = TRUE))
  overlap$strand <- strand(center)[overlap[,1]]
  overlap$name <- mcols(insertions)[overlap[,2],by]
  overlap <- transform(overlap, id=match(name, unique(name)))
  ids <- length(unique(overlap$name))

  #distance
  overlap$dist <- NA
  minus <- which(overlap$strand == "-")
  other <- which(overlap$strand != "-")
  overlap$dist[minus] <- start(center[overlap[minus,1]]) - start(insertions[overlap[minus,2]])
  overlap$dist[other] <- start(insertions[overlap[other,2]]) - start(center[overlap[other,1]])

  #Insertion Mat
  profile_mat <- tabulate2dCpp(x1 = overlap$id, y1 = overlap$dist, xmin = 1, xmax = ids, ymin = -flank, ymax = flank)
  colnames(profile_mat) <- unique(overlap$name)

  # insertion number fall into the flank of TSS
  profile <- rowSums(profile_mat)

  #normalize
  profile_mat_norm <- apply(profile_mat, 2,
                            function(x) x / max(mean(x[c(1:norm,(flank*2-norm+1):(flank*2+1))]), 0.5)) #Handles low depth cells
  profile_norm <- profile / mean(profile[c(1:norm,(flank*2-norm+1):(flank*2+1))])

  #smooth
  profile_mat_norm_smooth <- apply(profile_mat_norm, 2, function(x) zoo::rollmean(x, smooth, fill = 1))
  profile_norm_smooth <- zoo::rollmean(profile_norm, smooth, fill = 1)

  #enrichment
  max_finite <- function(x){
    suppressWarnings(max(x[is.finite(x)], na.rm=TRUE))
  }
  e_mat <- apply(profile_mat_norm_smooth, 2, function(x) max_finite(x[(flank-range):(flank+range)]))
  names(e_mat) <- colnames(profile_mat_norm_smooth)
  e <- max_finite(profile_norm_smooth[(flank-range):(flank+range)])

  #Summary
  df_mat <- data.frame(
    enrichment = e_mat,
    insertions = as.vector(table(mcols(insertions)[,by])[names(e_mat)]),
    insertionsWindow = as.vector(table(overlap$name)[names(e_mat)])
  )
  df_sum <- data.frame(bp = (-flank):flank, profile = profile, norm_profile = profile_norm, smooth_norm_profile = profile_norm_smooth, enrichment = e)
  rownames(df_sum) <-  NULL

  return(list(df = df_sum, dfall = df_mat, profileMat = profile_mat_norm, profileMatSmooth = profile_mat_norm_smooth))
}


