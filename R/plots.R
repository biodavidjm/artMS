# artMS PLOT FUNCTIONS
# 
# ------------------------------------------------------------------------------
#' @title Plot correlation distributions
#' 
#' @description Plot correlation distributions
#' @param MatrixCorrelations Matrix of correlations
#' @return A ggplot2 correlation plot
#' @keywords plot, correlation
#' artms_plotCorrelationDistribution()
#' @export
artms_plotCorrelationDistribution <- function(MatrixCorrelations){
  cor.data <- MatrixCorrelations[upper.tri(MatrixCorrelations, diag = FALSE)]  # we're only interested in one of the off-diagonals, otherwise there'd be duplicates
  cor.data <- as.data.frame(cor.data)  # that's how ggplot likes it
  colnames(cor.data) <- "pearson"
  
  g <- ggplot(data = cor.data, mapping = aes(x = pearson))
  g <- g + scale_x_continuous(breaks = seq(0, 1, by = 0.1)) 
  g <- g + geom_histogram(breaks=seq(0, 1, by =0.05), 
                          # col="red", 
                          aes(fill=..count..))
  g <- g + theme_minimal()
  g <- g + theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15)) 
  return(g)
}

# ------------------------------------------------------------------------------
#' @title Protein abundance dot plots
#' 
#' @description Protein abundance dot plots for each unique uniprot id
#' @param input_file The `-normalized.txt` output file from MSstats
#' @param output_file Wished output file name (add the `.pdf` extension)
#' @return A pdf file with each individual protein abundance plot for each
#' conditions
#' @keywords abundance, dotplots, plot
#' artms_dataPlots()
#' @export
artms_dataPlots <- function(input_file, output_file){
  
  data_mss = fread(input_file, integer64 = 'double')
  unique_subjects = unique(data_mss$PROTEIN)
  condition_length = length(unique(data_mss$GROUP_ORIGINAL))
  min_abu = min(data_mss$ABUNDANCE, na.rm = T)
  max_abu = max(data_mss$ABUNDANCE, na.rm=T)
  
  pdf(output_file, width = condition_length*1.5, height = 3)
  
  cat('PRINTING CONDITION PLOTS\n')
  for(subject in unique_subjects){
    subject_data = data_mss[PROTEIN==subject,]
    cat(sprintf('\t%s\n',subject))
    p = ggplot(data = subject_data, aes(x=SUBJECT_ORIGINAL,y=ABUNDANCE, colour=FEATURE))
    p = p + geom_point(size=2) + 
      facet_wrap(facets = ~ GROUP_ORIGINAL, drop = T, scales = 'free_x', ncol = condition_length) + 
      ylim(min_abu,max_abu) +
      theme(axis.text.x=element_text(angle=-90,hjust=1)) +
      guides(colour=FALSE) +
      xlab(NULL) +
      ggtitle(subject)
    print(p)
  }
  dev.off()
}

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
#' @param display Metric to be displayed (default: `log2fc`)
#' @return A heatmap in pdf format of the MSStats results using the 
#' log2 fold changes and the data used to generate the heatmap
#' @keywords heatmap, log2fc
#' artms_plotHeatmap()
#' @export
artms_plotHeatmap <- function(input_file, output_file, labels='*', cluster_cols=F, display='log2FC', lfc_lower=-2, lfc_upper=2, FDR=0.05){
  ## read input
  input <- read.delim(input_file, stringsAsFactors = F)
  
  ## select data points  by LFC & FDR criterium in single condition and adding corresponding data points from the other conditions
  sign_hits <- artms_significantHits(input,labels=labels,LFC=c(lfc_lower,lfc_upper),FDR=FDR)
  sign_labels <- unique(sign_hits$Label)
  cat(sprintf(">> SELECTED HITS FOR PLOTS WITH LFC BETWEEN %s AND %s AT %s FDR:\t%s\n",lfc_lower, lfc_upper, FDR, nrow(sign_hits)/length(sign_labels))) 
  
  ## REPRESENTING RESULTS AS HEATMAP
  ## plot heat map for all contrasts
  if(any(grepl('uniprot_genename',colnames(sign_hits)))){
    heat_labels <- paste(sign_hits$Protein,sign_hits$uniprot_genename,sep=' ')  
  }else{
    heat_labels <- sign_hits$Protein
  }
  
  heat_labels <- gsub('\\sNA$','',heat_labels)
  
  # Old PlotHeat function:
  heat_data = data.frame(sign_hits, heat_labels=heat_labels)
  
  ## create matrix from log2FC or p-value as user defined
  if(display=='log2FC'){
    # Issues with extreme_val later if we have Inf/-Inf values.
    if( sum(is.infinite(heat_data$log2FC)) > 0 ){
      idx <- is.infinite(heat_data$log2FC)
      heat_data$log2FC[ idx ] <- NA
    }
    heat_data_w = dcast(heat_labels ~ Label, data=heat_data, value.var='log2FC') 
  }else if(display=='adj.pvalue'){
    heat_data$adj.pvalue = -log10(heat_data$adj.pvalue+10^-16)  
    heat_data_w = dcast(heat_labels ~ Label, data=heat_data, value.var='adj.pvalue')  
  }else if(display=='pvalue'){
    heat_data$pvalue = -log10(heat_data$pvalue+10^-16)  
    heat_data_w = dcast(heat_labels ~ Label, data=heat_data, value.var='pvalue')  
  }
  
  ## try
  #gene_names = uniprot_to_gene_replace(uniprot_ac=heat_data_w$Protein)
  rownames(heat_data_w) = heat_data_w$heat_labels
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
  
  cat("--- Saving heatmap\n")
  pheatmap(heat_data_w, scale="none", cellheight=10, cellwidth=10, filename = output_file, color=colors_tot, breaks=seq(from=-extreme_val, to=extreme_val, by=bin_size), cluster_cols=cluster_cols, fontfamily="mono")
  return(heat_data_w)
}



# ------------------------------------------------------------------------------
#' @title Plot reproducibility plots
#' 
#' @description Plot reproducibility plots
#' @param data evidence data.frame
#' @return A reproducibility plot
#' @keywords plot, qc, quality, control
#' artms_plotReproducibilityEvidence()
#' @export
artms_plotReproducibilityEvidence <- function(data) {
  
  data <- data[c('Feature', 'Proteins', 'Intensity', 'Condition', 'BioReplicate', 'Run')]
  condi <- unique(data$Condition)
  
  # Progress bar
  pb <- txtProgressBar(min=0, max=length(condi), style=3)
  
  for(i in 1:length(condi)){
    eCondition <- condi[i]
    # Progress bar
    setTxtProgressBar(pb, i)
    
    # cat("\n\nCONDITION: ",eCondition,"\n##################\n")
    # cat("TECHNICAL REPLICAS\n---------------------------\n")
    # cat("- ", eCondition,"\n")
    
    conditionOne <- data[which(data$Condition == eCondition),]
    
    # FIRST CHECK FOR TECHNICAL REPLICAS
    bioreplicasAll <- unique(conditionOne$BioReplicate)
    
    for(eBioreplica in bioreplicasAll){
      
      # cat('\tChecking for technical replicas in ',eBioreplica, "\n")
      biorepli <- conditionOne[conditionOne$BioReplicate == eBioreplica,]
      here <- unique(biorepli$Run)
      
      if(length(here) > 1){ #Limit to 2 technical replicas: (length(here) > 1 & length(here) == 2)
        # We are expecting no more than 2 technical replicas. If there is more... it is worthy to double check
        # cat('\t\t>>Reproducibility for technical replicas of',eBioreplica,":")
        #Need to change the RUN number to letters (TR: TECHNICAL REPLICA)
        biorepli$TR <- biorepli$Run
        biorepli$TR[biorepli$TR == here[1]] <- 'tr1'
        biorepli$TR[biorepli$TR == here[2]] <- 'tr2'
        
        # Let's select unique features per TECHNICAL REPLICAS, and sum them up (the same feature might have been many differnt times)
        biorepliaggregated <- aggregate(Intensity~Feature+Proteins+Condition+BioReplicate+Run+TR, data = biorepli, FUN = sum)
        biorepliaggregated$Intensity <- log2(biorepliaggregated$Intensity)
        bdc <- dcast(data=biorepliaggregated, Feature+Proteins~TR, value.var = 'Intensity')
        
        # Get the number of proteins
        np <- dim(bdc)[1]
        corr_coef <- round(cor(bdc$tr1, bdc$tr2, use = "pairwise.complete.obs"), digits = 2)
        # cat("r:\t",corr_coef,"\n")
        p1 <- ggplot(bdc, aes(x=tr1, y = tr2))
        p1 <- p1 + geom_point() + geom_rug() + geom_density_2d(colour = 'lightgreen')
        p1 <- p1 + geom_smooth(colour = "green", fill = "lightblue", method = 'lm')
        p1 <- p1 + theme_light()
        p1 <- p1 + labs(title = paste("Reproducibility between Technical Replicas\nBioReplica:",eBioreplica, "  (n = ",np,", r = ",corr_coef,")"))
        print(p1)
      }else if(length(here) == 1){
        # cat("\t\tOnly one technical replica\n")
      }else{
        cat("\n>>MORE THAN TWO TECHNICAL REPLICAS IN THIS EXPERIMENTS? That is very strange.\n\n")
        stop("\nPlease, Check the keys files\n")
      }
    } # Checking the reproducibility between Technical Replicas
    
    
    # NOW BETWEEN BIOREPLICAS
    # Before comparing the different biological replicas, aggregate the technical replicas
    # cat("\nBIOLOGICAL REPLICAS\n---------------------------\n")
    # 
    # First choose the maximum for the technical replicas as before, but first check whether there are more than one
    if(length(here) > 1){
      conditionOne <- aggregate(Intensity~Feature+Proteins+Condition+BioReplicate+Run, data = conditionOne, FUN = sum)
      b <- aggregate(Intensity~Feature+Proteins+Condition+BioReplicate, data = conditionOne, FUN = mean)
    }else{
      b <- aggregate(Intensity~Feature+Proteins+Condition+BioReplicate, data = conditionOne, FUN = sum)
    }
    
    b$Intensity <- log2(b$Intensity)
    blist <- unique(b$BioReplicate)
    if(length(blist) > 1){ # We need at least TWO BIOLOGICAL REPLICAS
      
      to <- length(blist)-1
      
      for (i in 1:to) {
        
        j <- i+1
        
        for(k in j:length(blist)){
          
          br1 <- blist[i]
          br2 <- blist[k]
          
          # cat("\tChecking reproducibility between ",br1, "and ",br2 ,"\t")
          bcfinal <- dcast(data=b, Feature+Proteins~BioReplicate, value.var = 'Intensity')
          
          # Let's check the total number of peptides here...
          checkTotalNumber <- subset(bcfinal,select = c(br1, br2))
          # checkTotalNumber <- checkTotalNumber[complete.cases(checkTotalNumber),]
          
          npt <- dim(checkTotalNumber)[1]
          
          corr_coef <- round(cor(bcfinal[[br1]], bcfinal[[br2]], use = "pairwise.complete.obs"), digits = 2)
          # cat("r:\t",corr_coef,"\n")
          p2 <- ggplot(bcfinal, aes(x=bcfinal[[br1]], y = bcfinal[[br2]]))
          p2 <- p2 + geom_point()  + geom_rug() + geom_density_2d(colour = 'red')
          p2 <- p2 + geom_smooth(colour = "red", fill = "lightgreen", method = 'lm')
          p2 <- p2 + theme_light()
          p2 <- p2 + labs(title = paste("Peptide Reproducibility between Bioreplicas\n(condition:",eCondition,")\n",br1,"vs",br2,"(n =",npt," r = ",corr_coef,")"))
          p2 <- p2 + labs(x = br1)
          p2 <- p2 + labs(y = br2)
          print(p2)
        }
      }
      # cat("\n")
    }else{
      cat("\tONLY ONE BIOLOGICAL REPLICA AVAILABLE (plots are not possible)\n")
    }
  } # all the conditions
  # Close Progress bar
  close(pb)
}





