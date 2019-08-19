plot_GA_umap <- function(cds,
                         cicero_gene_activities,
                         gene,
                         dotSize = 1,
                         log = TRUE){
  umap_df <- as.data.frame(colData(cds))
  gene_pos <- which(row.names(cicero_gene_activities) %in% gene)

  umap_df$score <- cicero_gene_activities[gene_pos,]
  if (log){
    umap_df$score <- log(umap_df$score + 1)
  }

  p <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
    geom_point(aes(color=score), size = dotSize) + scale_color_viridis(option="magma") +
    theme_bw()
  return(p)
}
