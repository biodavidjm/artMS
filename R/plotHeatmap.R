# ------------------------------------------------------------------------------
#' @title Outputs a heatmap of the MSStats results created using the log2fold 
#' changes
#' 
#' @description Outputs a heatmap of the MSStats results created using the 
#' log2 fold changes.
#' @param input_file MSstats `results.txt` file and location
#' @param output_file Output file name (pdf format) and location
#' @param labels Vector of uniprot ids if only specific labes would like to
#' be plotted. Default: all labels
#' @param cluster_cols `True` or `False` to cluster columns. Default: FALSE
#' @param lfc_lower Lower limit for the log2fc. Default: -2
#' @param lfc_upper Upper limit for the log2fc. Default: +2
#' @param FDR Upper limit false discovery rate. Default: 0.05
#' @return A heatmap in pdf format of the MSStats results using the 
#' log2 fold changes.
#' @keywords heatmap, log2fc
#' artms_plotHeatmap()
#' @export
artms_plotHeatmap <- function(input_file, output_file, labels='*', cluster_cols=F, display='log2FC', lfc_lower=-2, lfc_upper=2, FDR=0.05){
  ## read input
  input <- read.delim(input_file, stringsAsFactors = F)
  
  ## select data points  by LFC & FDR criterium in single condition and adding corresponding data points from the other conditions
  sign_hits <- artms_significantHits(input,labels=labels,LFC=c(lfc_lower,lfc_upper),FDR=FDR)
  sign_labels <- unique(sign_hits$Label)
  cat(sprintf("\tSELECTED HITS FOR PLOTS WITH LFC BETWEEN %s AND %s AT %s FDR:\t%s\n",lfc_lower, lfc_upper, FDR, nrow(sign_hits)/length(sign_labels))) 
  
  ## REPRESENTING RESULTS AS HEATMAP
  ## plot heat map for all contrasts
  if(any(grepl('uniprot_genename',colnames(sign_hits)))){
    heat_labels <- paste(sign_hits$Protein,sign_hits$uniprot_genename,sep=' ')  
  }else{
    heat_labels <- sign_hits$Protein
  }
  
  heat_labels <- gsub('\\sNA$','',heat_labels)
  heat_data_w <- plotHeat(mss_F=sign_hits, out_file=output_file, names=heat_labels, cluster_cols=cluster_cols, display=display)  
}


# ------------------------------------------------------------------------------
#' @title Heatmap of significant values
#' 
#' @description heatmap plot to represent proteins with significant changes
#' @param mss_F data.frame with the significant values (log2fc, pvalues)
#' @param out_file Name for the output
#' @param labelOrder Vector with the particular order for the IDs (default, 
#' `NULL` no order)
#' @param names Type of ID used. Default is `Protein`` (uniprot entry id). 
#' Soon will be possible to use 'Gene' name ids.
#' @param cluster_cols Select whether to cluster the columns. Options: `T` 
#' or `F`. Default `T`.
#' @param display Value used to genarate the heatmaps. Options: `log2FC`, 
#' `adj.pvalue`, `pvalue`. Default: `log2FC`
#' @return A heatmap of significant values
#' @keywords significant, heatmap
#' plotHeat()
#' @export
plotHeat <- function(mss_F, out_file, labelOrder=NULL, names='Protein', cluster_cols=F, display='log2FC'){
  heat_data = data.frame(mss_F, names=names)
  
  ## create matrix from log2FC or p-value as user defined
  if(display=='log2FC'){
    # Issues with extreme_val later if we have Inf/-Inf values.
    if( sum(is.infinite(heat_data$log2FC)) > 0 ){
      idx <- is.infinite(heat_data$log2FC)
      heat_data$log2FC[ idx ] <- NA
    }
    heat_data_w = dcast(names ~ Label, data=heat_data, value.var='log2FC') 
  }else if(display=='adj.pvalue'){
    heat_data$adj.pvalue = -log10(heat_data$adj.pvalue+10^-16)  
    heat_data_w = dcast(names ~ Label, data=heat_data, value.var='adj.pvalue')  
  }else if(display=='pvalue'){
    heat_data$pvalue = -log10(heat_data$pvalue+10^-16)  
    heat_data_w = dcast(names ~ Label, data=heat_data, value.var='pvalue')  
  }
  
  ## try
  #gene_names = uniprot_to_gene_replace(uniprot_ac=heat_data_w$Protein)
  rownames(heat_data_w) = heat_data_w$names
  heat_data_w = heat_data_w[,-1]
  heat_data_w[is.na(heat_data_w)]=0
  max_val = ceiling(max(heat_data_w))
  min_val = floor(min(heat_data_w))
  extreme_val = max(max_val, abs(min_val))
  if(extreme_val %% 2 != 0) extreme_val=extreme_val+1
  bin_size=2
  signed_bins = (extreme_val/bin_size)
  colors_neg = rev(colorRampPalette(RColorBrewer::brewer.pal("Blues",n=extreme_val/bin_size))(signed_bins))
  colors_pos = colorRampPalette(RColorBrewer::brewer.pal("Reds",n=extreme_val/bin_size))(signed_bins)
  colors_tot = c(colors_neg, colors_pos)
  
  if(is.null(labelOrder)){
    cat("\t Saving heatmap\n")
    pheatmap(heat_data_w, scale="none", cellheight=10, cellwidth=10, filename =out_file, color=colors_tot, breaks=seq(from=-extreme_val, to=extreme_val, by=bin_size), cluster_cols=cluster_cols, fontfamily="mono")  
  }else{
    heat_data_w = heat_data_w[,labelOrder]
    pheatmap(heat_data_w, scale="none", cellheight=10, cellwidth=10, filename=out_file, color=colors_tot, breaks=seq(from=-extreme_val, to=extreme_val, by=bin_size), cluster_cols=cluster_cols, fontfamily="mono")
  }
  return(heat_data_w)
}