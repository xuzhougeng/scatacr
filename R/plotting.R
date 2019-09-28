#' Plot TSS Profile
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_hex
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom viridis scale_fill_viridis
#' @importFrom stats complete.cases
#'
#' @param tssSingles object return from getTssProfile
#' @param filterFrags threshold of fragments of cell
#' @param filterTSS enrichment score of TSS
#' @param outFile graphics file path
#' @rdname plotTssProfile
#' @export
#'
plotTssProfile <- function(tssSingles,
                           filterFrags = 1000,
                           filterTSS = 8,
                           outFile = NULL){

  tssSingles$cellCall <- 0
  tssSingles$cellCall[tssSingles$uniqueFrags >= filterFrags &
                        tssSingles$enrichment >= filterTSS] <- 1


  tssSingles <- tssSingles[complete.cases(tssSingles),]
  nPass  <- sum(tssSingles$cellCall==1)
  nTotal <- sum(tssSingles$uniqueFrags >= filterFrags)

  p <- ggplot(tssSingles[tssSingles$uniqueFrags > 500,],
              aes(x = log10(uniqueFrags), y = enrichment)) +
    geom_hex(bins = 100) +
    theme_bw() + scale_fill_viridis() +
    xlab("log10 Unique Fragments") +
    ylab("TSS Enrichment") +
    geom_hline(yintercept = filterTSS, lty = "dashed") +
    geom_vline(xintercept = log10(filterFrags), lty = "dashed") +
    ggtitle(sprintf("Pass Rate : %s of %s (%s)", nPass, nTotal, round(100*nPass/nTotal,2)))

  if( is.null(outFile)){
    return(p)
  } else{
    pdf(outFile)
    print(p)
    dev.off()
    return(p)
  }

}


#' plot gene activity score in UMAP
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme_bw
#' @importFrom viridis scale_color_viridis
#'
#' @rdname plot_GA_umap
#' @export
plot_GA_umap <- function(cds,
                         cicero_gene_activities,
                         gene,
                         dotSize = 1,
                         scale = 1000000){
  umap_df <- as.data.frame(colData(cds))
  gene_pos <- which(row.names(cicero_gene_activities) %in% gene)

  umap_df$score <- cicero_gene_activities[gene_pos,]

  umap_df$score <- log(umap_df$score * scale + 1)


  p <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
    geom_point(aes(color=score), size = dotSize) +
    scale_color_viridis(option="magma") +
    theme_bw()
  return(p)
}
