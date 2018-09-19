# ------------------------------------------------------------------------------
#' @title Volcano plot (log2fc / pvalues)
#' 
#' @description It generates a scatter-plot used to quickly identify changes
#' @param mss_results_sel (data.frame) Selected MSstats results
#' @param lfc_upper (numeric) log2fc upper threshold (positive value)
#' @param lfc_lower (numeric) log2fc lower threshold (negative value)
#' @param FDR (numeric) False Discovery Rate threshold
#' @param file_name (char) Name for the output file
#' @param PDF (logical) Option to generate pdf format. Default: `T`
#' @param decimal_threshold (numeric) Decimal threshold for the pvalue. 
#' Default: 16 (10^-16)
#' @keywords plot, volcano
#' @return (pdf) of a volcano plot
#' @examples \donttest{
#' # Open the msstats results
#' mss <- read.delim("resultsQuant/a549-PB1-results.txt", stringsAsFactors = FALSE)
#' # And generate volcano plot
#' artms_volcanoPlot(mss_results_sel = mss,
#'                   lfc_upper = 1, 
#'                   lfc_lower = -1, 
#'                   FDR = 0.05, 
#'                   file_name = "a549-PB1-results-volcanoPlot.pdf", 
#'                   PDF = TRUE)
#' }
#' @export
artms_volcanoPlot <- function(mss_results_sel, 
                              lfc_upper, 
                              lfc_lower, 
                              FDR, 
                              file_name='', 
                              PDF= TRUE, 
                              decimal_threshold=16){
  
  # handle cases where log2FC is Inf. There are no pvalues or other information for these cases :(
  # Issues with extreme_val later if we have Inf/-Inf values.
  if( sum(is.infinite(mss_results_sel$log2FC)) > 0 ){
    idx <- is.infinite(mss_results_sel$log2FC)
    mss_results_sel$log2FC[ idx ] <- NA
  }
  
  min_x = -ceiling(max(abs(mss_results_sel$log2FC), na.rm= TRUE))
  max_x = ceiling(max(abs(mss_results_sel$log2FC), na.rm= TRUE))
  # Deal with special cases in the data where we have pvalues = Inf,NA,0
  if( sum(is.na(mss_results_sel$adj.pvalue))>0 ) mss_results_sel <- mss_results_sel[!is.na(mss_results_sel$adj.pvalue),]
  if(nrow(mss_results_sel[mss_results_sel$adj.pvalue == 0 | mss_results_sel$adj.pvalue == -Inf,]) > 0) mss_results_sel[!is.na(mss_results_sel$adj.pvalue) & (mss_results_sel$adj.pvalue == 0 | mss_results_sel$adj.pvalue == -Inf),]$adj.pvalue = 10^-decimal_threshold
  max_y = ceiling(-log10(min(mss_results_sel[mss_results_sel$adj.pvalue > 0,]$adj.pvalue, na.rm= TRUE))) + 1
  
  l = length(unique(mss_results_sel$Label))
  w_base = 7
  h_base = 7
  
  if(l<=2){
    w=w_base*l 
  }else{
    w=w_base*2
  }
  h = h_base*ceiling(l/2)
  
  if(PDF) pdf(file_name, width=w, height=h)
  p = ggplot(mss_results_sel, aes(x=log2FC,y=-log10(adj.pvalue)))
  print(p + geom_point(colour='grey') + 
          geom_point(data = mss_results_sel[mss_results_sel$adj.pvalue <= FDR & mss_results_sel$log2FC>=lfc_upper,], aes(x=log2FC,y=-log10(adj.pvalue)), colour='red', size=2) +
          geom_point(data = mss_results_sel[mss_results_sel$adj.pvalue <= FDR & mss_results_sel$log2FC<=lfc_lower,], aes(x=log2FC,y=-log10(adj.pvalue)), colour='blue', size=2) +
          geom_vline(xintercept=c(lfc_lower,lfc_upper), lty='dashed') + 
          geom_hline(yintercept=-log10(FDR), lty='dashed') + 
          xlim(min_x,max_x) + 
          ylim(0,max_y) + 
          facet_wrap(facets = ~Label, ncol = 2, scales = 'fixed')) 
  if(PDF) dev.off()
}