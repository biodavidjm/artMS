#' @title Create Volcano Plots from MSstats results
#' @description This function creates volcanon plots of the MSstats results, allowing the user to define regions of interest (confident hits) based on the Log2FC and adjusted p-value. A ggplot object will be returned. 
#' @param mss_results_sel The filepath to the MSstats results file (txt tab delimited file).
#' @param file_name The filepath to the intended output file (txt tab delimited file).
#' @param lfc_upper The starting point of the upper region on interest based on the Log2FC (default = 2).
#' @param lfc_lower The starting point of the lower region on interest based on the Log2FC (default = -2).
#' @param FDR The adjusted p-value cutoff to define the region of interest (default = 0.05).
#' @param PDF Whether or not to write ou a pdf of the results or not. (default = F).
#' @param decimal_threshold Apply a decimal threshold. (default = 16).
#' @keywords volcano, MSstats
#' volcanoPlot()
volcanoPlot = function(mss_results_sel, lfc_upper = 2, lfc_lower = -2, FDR = 0.05, file_name = "", PDF = F, decimal_threshold = 16) {
    
    mss_results_sel = checkIfFile(mss_results_sel, is.evidence = FALSE)
    
    # create a column of useful text to pop up on mouseover when used with ggplotly
    mss_results_sel$Text <- paste0(mss_results_sel$Gene.names, "</br>", mss_results_sel$Protein)
    
    # handle cases where log2FC is Inf. There are no pvalues or other information for these cases :( Issues with extreme_val later
    # if we have Inf/-Inf values.
    if (sum(is.infinite(mss_results_sel$log2FC)) > 0) {
        idx <- is.infinite(mss_results_sel$log2FC)
        mss_results_sel$log2FC[idx] <- NA
    }
    
    min_x = -ceiling(max(abs(mss_results_sel$log2FC), na.rm = T))
    max_x = ceiling(max(abs(mss_results_sel$log2FC), na.rm = T))
    # Deal with special cases in the data where we have pvalues = Inf,NA,0
    if (sum(is.na(mss_results_sel$adj.pvalue)) > 0) 
        mss_results_sel <- mss_results_sel[!is.na(mss_results_sel$adj.pvalue), ]
    if (nrow(mss_results_sel[mss_results_sel$adj.pvalue == 0 | mss_results_sel$adj.pvalue == -Inf, ]) > 0) 
        mss_results_sel[!is.na(mss_results_sel$adj.pvalue) & (mss_results_sel$adj.pvalue == 0 | mss_results_sel$adj.pvalue == -Inf), 
            ]$adj.pvalue = 10^-decimal_threshold
    max_y = ceiling(-log10(min(mss_results_sel[mss_results_sel$adj.pvalue > 0, ]$adj.pvalue, na.rm = T))) + 1
    
    l = length(unique(mss_results_sel$Label))
    w_base = 7
    h_base = 7
    
    if (l <= 2) {
        w = w_base * l
    } else {
        w = w_base * 2
    }
    h = h_base * ceiling(l/2)
    
    if (PDF) 
        pdf(file_name, width = w, height = h)
    grey_idx <- (mss_results_sel$adj.pvalue > FDR) & (mss_results_sel$log2FC < lfc_upper) & (mss_results_sel$log2FC > lfc_lower)
    p = ggplot(mss_results_sel[grey_idx, ], aes(x = log2FC, y = -log10(adj.pvalue), text = Text))
    p = p + geom_point(colour = "grey") + geom_point(data = mss_results_sel[mss_results_sel$adj.pvalue <= FDR & mss_results_sel$log2FC >= 
        lfc_upper, ], aes(x = log2FC, y = -log10(adj.pvalue)), colour = "red", size = 2) + geom_point(data = mss_results_sel[mss_results_sel$adj.pvalue <= 
        FDR & mss_results_sel$log2FC <= lfc_lower, ], aes(x = log2FC, y = -log10(adj.pvalue)), colour = "blue", size = 2) + geom_vline(xintercept = c(lfc_lower, 
        lfc_upper), lty = "dashed") + geom_hline(yintercept = -log10(FDR), lty = "dashed") + xlim(min_x, max_x) + ylim(0, max_y) + 
        facet_wrap(facets = ~Label, ncol = 2, scales = "fixed")
    print(p)
    if (PDF) 
        dev.off()
    return(p)
}





