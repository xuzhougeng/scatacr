#' plot gene activity score in UMAP
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme_bw
#' @importFrom virids scale_color_viridis
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
