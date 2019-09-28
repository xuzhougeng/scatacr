
#' Generate insertion count matrixy by feature
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors DataFrame
#' @importFrom S4Vectors mcols
#'
#' @param feature GenomicRanges object  to store feature information
#' @param fragments GenomicRanges object
#' @param by metadata in mcols
#' @param ... Arguments passed to other methods
#'
#' @return An list with sparseMatrix, FRIP and total fragments
#'
#' @rdname countInsertions
#' @export
countInsertions <- function(fragments, feature,  by = "RG", ...){
  #Count By Fragments Insertions
  inserts <- c(
    GRanges(seqnames = seqnames(fragments),
            ranges = IRanges(start(fragments), start(fragments)),
            RG = mcols(fragments)[,by]),
    GRanges(seqnames = seqnames(fragments),
            ranges = IRanges(end(fragments), end(fragments)),
            RG = mcols(fragments)[,by])
  )
  by <- "RG"
  overlapDF <- DataFrame(findOverlaps(feature, inserts,
                                      ignore.strand = TRUE,
                                      maxgap=-1L, minoverlap=0L, type = "any"))
  overlapDF$name <- mcols(inserts)[overlapDF[, 2], by]
  overlapTDF <- transform(overlapDF, id = match(name, unique(name)))
  #Calculate Overlap Stats
  inPeaks <- table(overlapDF$name)
  total <- table(mcols(inserts)[, by])
  total <- total[names(inPeaks)]
  frip <- inPeaks / total
  #Summarize
  sparseM <- Matrix::sparseMatrix(
    i = overlapTDF[, 1],
    j = overlapTDF[, 4],
    x = rep(1, nrow(overlapTDF)),
    dims = c(length(feature), length(unique(overlapDF$name))))
  colnames(sparseM) <- unique(overlapDF$name)
  total <- total[colnames(sparseM)]
  frip <- frip[colnames(sparseM)]
  out <- list(counts = sparseM, frip = frip, total = total)
  return(out)
}
