#' @import ggplot2
#' @import ggdendro
#' @import plotly
#' @import reshape2

# library(ggplot2)
# library(ggdendro)
# library(plotly)
# library(reshape2)
results_file <- '~/Box Sync/projects/krogan/rkaake/apobec/msstats/results/20170227_equalizedMedians/20170227-results.txt'
pval <- 0.05

#' @title Create Heatmap from MSStat Results
#' @description This will read in the MSStats results file (original long form) and output a heatmap that the user can interact with by hovering the mouse over it.
#' @param results_file The filepath to the MSstats results file in the original long format.
#' @keywords MSStats
#' resultsHeatmap()
#' @export
resultsHeatmap <- function(results_file, pval=0.05, prot_file=FALSE){
  results <- checkIfFile(results_file, is.evidence=FALSE)

  if(prot_file){
    ##### FILL THIS IN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  }else{
    #remove non significant cases
    cat('NO PROTEIN FILE PROVIDED TO FILTER PROTEINS WITH. USING WHOLE DATA SET')
    # get a list of protiens that are significant and meet the log2fc criterion. Also log2fc must not = +/-Inf
    prots <- unique(results$Protein[which( (results$adj.pvalue < pval) & !is.infinite(results$log2FC) & !is.na(results$log2FC) )])

  }

  # pull out all cases of these proteins in the data set
  results.thresh <- results[which(results$Protein %in% prots),]
  
  # fill in +/-Inf values with the Max/Min value
  x.max <- max(results.thresh$log2FC[which(!is.infinite(results.thresh$log2FC))], na.rm = T)
  x.min <- min(results.thresh$log2FC[which(!is.infinite(results.thresh$log2FC))], na.rm = T)
  results.thresh$log2FC[which(results.thresh$log2FC == Inf)] <- x.max
  results.thresh$log2FC[which(results.thresh$log2FC == -Inf)] <- x.min
  
  
  # convert to wide format
  x <- results.thresh
  x <- dcast(data=x, Protein+Gene.names~Label, value.var='log2FC', fill=0)
  row.names(x) = paste( x$Protein, x$Gene.names, sep='-')
  x$Protein = x$Gene.names = c()
  x.names = row.names(x)

  #dendogram data
  # x <- as.matrix(scale(x))
  x <- as.matrix(x)
  row.names(x) = x.names
  dd.col <- as.dendrogram(hclust(dist(x)))
  dd.row <- as.dendrogram(hclust(dist(t(x))))
  dx <- dendro_data(dd.row)
  dy <- dendro_data(dd.col)

  # Change order of heatmap to reflect dendrogram clustering
  col.ord <- order.dendrogram(dd.col)
  row.ord <- order.dendrogram(dd.row)
  xx <- scale(x)[col.ord, row.ord]
  # change names
  xx_names <- attr(xx, "dimnames")
  df <- as.data.frame(xx)
  colnames(df) <- xx_names[[2]]
  
  df$Proteins <- xx_names[[1]]

  df$Proteins <- with(df, factor(Proteins, levels=Proteins, ordered=TRUE))

  # put data into ggplot format
  mdf <- reshape2::melt(df, id.vars="Proteins")
  p <- ggplot(mdf, aes(x = variable, y = Proteins)) + geom_tile(aes(fill = value), colour = "white") + theme(axis.text.x = element_text(angle=45, hjust=1)) + labs(title="Log2 Fold Changes") + xlab("Protein") + ylab("Contrast") +
    scale_fill_gradient(low="blue",high="red")
  ggplotly(p)

}


































