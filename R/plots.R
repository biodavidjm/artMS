# artMS PLOT FUNCTIONS
#
# ------------------------------------------------------------------------------
#' @title Plot correlation distributions
#'
#' @description Plot correlation distributions
#' @param MatrixCorrelations (matrix) of correlations
#' @return A ggplot2 correlation plot
#' @keywords plot, correlation
.artms_plotCorrelationDistribution <- function(MatrixCorrelations) {
  # we're only interested in one of the off-diagonals,
  # otherwise there'd be duplicates
  cor.data <-
    MatrixCorrelations[upper.tri(MatrixCorrelations, diag = FALSE)]
  cor.data <- as.data.frame(cor.data)
  colnames(cor.data) <- "pearson"
  
  g <- ggplot(data = cor.data, mapping = aes(x = pearson))
  g <- g + scale_x_continuous(breaks = seq(0, 1, by = 0.1))
  g <- g + geom_histogram(breaks = seq(0, 1, by = 0.05),
                          # col="red",
                          aes(fill = ..count..))
  g <- g + theme_minimal()
  g <-
    g + theme(axis.text.x = element_text(size = 15),
              axis.text.y = element_text(size = 15))
  return(g)
}

# ------------------------------------------------------------------------------
#' @title Individual Normalized abundance dot plots for every protein
#'
#' @description Protein abundance dot plots for each unique uniprot id. It can
#' take a long time
#' @param input_file (char) File path and name to the `-normalized.txt` output
#' file from MSstats
#' @param output_file (char) Output file (path) name (add the `.pdf` extension)
#' @return (pdf) file with each individual protein abundance plot for each
#' conditions
#' @keywords abundance, dotplots, plot
#' @examples \donttest{
#' artms_dataPlots(input_file = "results/ab-results-mss-normalized.txt",
#'                output_file = "results/ab-results-mss-normalized.pdf")
#' }
#' @export
artms_dataPlots <- function(input_file, output_file) {
  data_mss = fread(input_file, integer64 = 'double')
  unique_subjects <- unique(data_mss$PROTEIN)
  condition_length <- length(unique(data_mss$GROUP_ORIGINAL))
  min_abu <- min(data_mss$ABUNDANCE, na.rm = TRUE)
  max_abu <- max(data_mss$ABUNDANCE, na.rm = TRUE)
  
  pdf(output_file, width = condition_length * 1.5, height = 3)
  cat('>> PRINTING CONDITION PLOTS for every protein\n')
  for (subject in unique_subjects) {
    subject_data <- data_mss[PROTEIN == subject, ]
    cat(sprintf('%s ', subject))
    p <-
      ggplot(data = subject_data,
             aes(x = SUBJECT_ORIGINAL, y = ABUNDANCE, colour = FEATURE))
    p <- p + geom_point(size = 2) +
      facet_wrap(
        facets = ~ GROUP_ORIGINAL,
        drop = TRUE,
        scales = 'free_x',
        ncol = condition_length
      ) +
      ylim(min_abu, max_abu) +
      theme(axis.text.x = element_text(angle = -90, hjust = 1)) +
      guides(colour = FALSE) +
      xlab(NULL) +
      ggtitle(subject)
    print(p)
  }
  cat("--- Done!\n")
  garbarge <- dev.off()
}

# ------------------------------------------------------------------------------
#' @title Heatmap of significant values
#'
#' @description heatmap plot to represent proteins with significant changes
#' @param mss_F (data.frame) with the significant values (log2fc, pvalues)
#' @param out_file (char) Name for the output
#' @param labelOrder (vector) Vector with the particular order for the IDs
#' (default, `NULL` no order)
#' @param names (char) Type of ID used. Default is `Protein` (uniprot entry id).
#' Soon will be possible to use 'Gene' name ids.
#' @param cluster_cols (logical) Select whether to cluster the columns.
#' Options: `T` or `F`. Default `T`.
#' @param display (char) Value used to genarate the heatmaps. Options:
#' - `log2FC` (default)
#' - `adj.pvalue`
#' - `pvalue`
#' @return A heatmap of significant values
#' @keywords significant, heatmap
.artms_plotHeat <-
  function(mss_F,
           out_file,
           labelOrder = NULL,
           names = 'Protein',
           cluster_cols = FALSE,
           display = 'log2FC') {
    heat_data = data.frame(mss_F, names = names)
    
    ## create matrix from log2FC or p-value as user defined
    if (display == 'log2FC') {
      # Issues with extreme_val later if we have Inf/-Inf values.
      if (sum(is.infinite(heat_data$log2FC)) > 0) {
        idx <- is.infinite(heat_data$log2FC)
        heat_data$log2FC[idx] <- NA
      }
      heat_data_w = dcast(names ~ Label, data = heat_data, value.var = 'log2FC')
    } else if (display == 'adj.pvalue') {
      heat_data$adj.pvalue = -log10(heat_data$adj.pvalue + 10 ^ -16)
      heat_data_w = dcast(names ~ Label, 
                          data = heat_data, value.var = 'adj.pvalue')
    } else if (display == 'pvalue') {
      heat_data$pvalue = -log10(heat_data$pvalue + 10 ^ -16)
      heat_data_w = dcast(names ~ Label, data = heat_data, value.var = 'pvalue')
    }
    
    ## try
    #gene_names = uniprot_to_gene_replace(uniprot_ac=heat_data_w$Protein)
    rownames(heat_data_w) = heat_data_w$names
    heat_data_w = heat_data_w[, -1]
    heat_data_w[is.na(heat_data_w)] = 0
    max_val = ceiling(max(heat_data_w))
    min_val = floor(min(heat_data_w))
    extreme_val = max(max_val, abs(min_val))
    if (extreme_val %% 2 != 0)
      extreme_val = extreme_val + 1
    bin_size = 2
    signed_bins = (extreme_val / bin_size)
    colors_neg = rev(colorRampPalette(brewer.pal("Blues", n = extreme_val /
                                                   bin_size))(signed_bins))
    colors_pos = colorRampPalette(
      brewer.pal("Reds", n = extreme_val / bin_size))(signed_bins)
    colors_tot = c(colors_neg, colors_pos)
    
    if (is.null(labelOrder)) {
      pheatmap(
        heat_data_w,
        scale = "none",
        cellheight = 10,
        cellwidth = 10,
        filename = out_file,
        color = colors_tot,
        breaks = seq(
          from = -extreme_val,
          to = extreme_val,
          by = bin_size
        ),
        cluster_cols = cluster_cols,
        fontfamily = "mono"
      )
      cat("--- Heatmap is out\n")
    } else{
      heat_data_w <- heat_data_w[, labelOrder]
      pheatmap(
        heat_data_w,
        scale = "none",
        cellheight = 10,
        cellwidth = 10,
        filename = out_file,
        color = colors_tot,
        breaks = seq(
          from = -extreme_val,
          to = extreme_val,
          by = bin_size
        ),
        cluster_cols = cluster_cols,
        fontfamily = "mono"
      )
      cat("--- Heatmap is out\n")
    }
    
    return(heat_data_w)
  }

# ------------------------------------------------------------------------------
#' @title Outputs a heatmap of the MSStats results created using the log2fold
#' changes
#'
#' @description Heatmap of the Relative Quantifications (MSStats results)
#' @param input_file (char) MSstats `results.txt` file and location (or
#' data.frame of resuts)
#' @param output_file (char) Output file name (pdf format) and location.
#' Default:"quantifications_heatmap.pdf"
#' @param specie (char). Specie name to be able to add the Gene name. To find
#' out more about the supported species check `?artms_mapUniprot2entrezGeneName`
#' @param labels (vector) of uniprot ids if only specific labes would like to
#' be plotted. Default: all labels
#' @param cluster_cols (boolean) `True` or `False` to cluster columns.
#' Default: FALSE
#' @param lfc_lower (int) Lower limit for the log2fc. Default: -2
#' @param lfc_upper (int) Upper limit for the log2fc. Default: +2
#' @param whatPvalue (char) `pvalue` or `adj.pvalue` (default)
#' @param FDR (int) Upper limit false discovery rate (or pvalue). Default: 0.05
#' @param display Metric to be displayed. Options:
#' - `log2fc` (default)
#' - `adj.pvalue`
#' - `pvalue`
#' @return (pdf or ggplot2 object) heatmap of the MSStats results using the
#' selected metric
#' @keywords heatmap, log2fc
#' artms_plotHeatmapQuant(input_file = artms_data_ph_msstats_results,
#'                        specie = "human",
#'                        output_file = NULL,
#'                        whatPvalue = "pvalue",
#'                        lfc_lower = -1,
#'                        lfc_upper = 1)
#' @export
artms_plotHeatmapQuant <- function(input_file,
                                   output_file = "quantifications_heatmap.pdf",
                                   specie,
                                   labels = '*',
                                   cluster_cols = FALSE,
                                   display = 'log2FC',
                                   lfc_lower = -2,
                                   lfc_upper = 2,
                                   whatPvalue = "adj.pvalue",
                                   FDR = 0.05) {
  input <- .artms_checkIfFile(input_file)
  
  ## select data points  by LFC & FDR criterium in single condition and
  ## adding corresponding data points from the other conditions
  sign_hits <- .artms_significantHits(
    input,
    labels = labels,
    LFC = c(lfc_lower, lfc_upper),
    whatPvalue = whatPvalue,
    FDR = FDR
  )
  sign_hits <- sign_hits[complete.cases(sign_hits$log2FC), ]
  sign_hits <- sign_hits[is.finite(sign_hits$log2FC), ]
  if (dim(sign_hits)[1] == 0) {
    stop("--- NOT ENOUGH SIGNIFICANT HITS!")
  }
  
  sign_labels <- unique(sign_hits$Label)
  cat(
    sprintf(
      ">> TOTAL NUMBER OF SELECTED HITS FOR PLOTS WITH LFC BETWEEN %s AND %s AT %s FDR:%s\n",
      lfc_lower,
      lfc_upper,
      FDR,
      nrow(sign_hits) / length(sign_labels)
    )
  )
  
  suppressMessages(
    sign_hits <- artms_annotationUniprot(
      data = sign_hits,
      columnid = "Protein",
      sps = specie
    )
  )
  
  ## REPRESENTING RESULTS AS HEATMAP
  ## plot heat map for all contrasts
  if (any(grepl('Gene', colnames(sign_hits)))) {
    heat_labels <- paste(sign_hits$Protein, sign_hits$Gene, sep = '_')
  } else{
    heat_labels <- sign_hits$Protein
  }
  
  heat_labels <- gsub('\\sNA$', '', heat_labels)
  
  # Old PlotHeat function:
  heat_data <- data.frame(sign_hits, heat_labels = heat_labels)
  
  ## create matrix from log2FC or p-value as user defined
  if (display == 'log2FC') {
    # Issues with extreme_val later if we have Inf/-Inf values.
    if (sum(is.infinite(heat_data$log2FC)) > 0) {
      idx <- is.infinite(heat_data$log2FC)
      heat_data$log2FC[idx] <- NA
    }
    heat_data_w <-
      data.table::dcast(heat_labels ~ Label, 
                        data = heat_data, 
                        value.var = 'log2FC')
  } else if (display == 'adj.pvalue') {
    heat_data$adj.pvalue <- -log10(heat_data$adj.pvalue + 10 ^ -16)
    heat_data_w <-
      data.table::dcast(heat_labels ~ Label, 
                        data = heat_data, 
                        value.var = 'adj.pvalue')
  } else if (display == 'pvalue') {
    heat_data$pvalue <- -log10(heat_data$pvalue + 10 ^ -16)
    heat_data_w <-
      data.table::dcast(heat_labels ~ Label, 
                        data = heat_data, 
                        value.var = 'pvalue')
  }
  
  #gene_names <- uniprot_to_gene_replace(uniprot_ac=heat_data_w$Protein)
  rownames(heat_data_w) <- heat_data_w$heat_labels
  heat_data_w <- heat_data_w[, -1]
  heat_data_w[is.na(heat_data_w)] = 0
  max_val <- ceiling(max(heat_data_w))
  min_val <- floor(min(heat_data_w))
  extreme_val <- max(max_val, abs(min_val))
  if (extreme_val %% 2 != 0)
    extreme_val = extreme_val + 1
  bin_size = 2
  signed_bins <- (extreme_val / bin_size)
  colors_neg <-
    rev(colorRampPalette(RColorBrewer::brewer.pal("Blues", n = extreme_val /
                                                    bin_size))(signed_bins))
  colors_pos <-
    colorRampPalette(RColorBrewer::brewer.pal("Reds", n = extreme_val / bin_size))(signed_bins)
  colors_tot <- c(colors_neg, colors_pos)
  
  if (!is.null(output_file)) {
    pheatmap(
      heat_data_w,
      scale = "none",
      cellheight = 10,
      cellwidth = 10,
      filename = output_file,
      color = colors_tot,
      breaks = seq(
        from = -extreme_val,
        to = extreme_val,
        by = bin_size
      ),
      cluster_cols = cluster_cols,
      fontfamily = "mono"
    )
    cat("--- Heatmap done\n")
  } else{
    pheatmap(
      heat_data_w,
      scale = "none",
      cellheight = 10,
      cellwidth = 10,
      color = colors_tot,
      breaks = seq(
        from = -extreme_val,
        to = extreme_val,
        by = bin_size
      ),
      cluster_cols = cluster_cols,
      fontfamily = "mono"
    )
  }
}

# ------------------------------------------------------------------------------
#' @title Generate reproducibility plots based on raw intentities
#' (from the evidence file)
#'
#' @description Generate reproducibility plots based on raw intentities
#' (log tranformed) from the evidence file
#' @param data (data.frame) clean processed evidence
#' @return (pdf) A reproducibility plot based on evidence values
#' @keywords internal, plot, qc, quality, control
.artms_plotReproducibilityEvidence <- function(data) {
  data <-
    data[c('Feature',
           'Proteins',
           'Intensity',
           'Condition',
           'BioReplicate',
           'Run')]
  condi <- unique(data$Condition)
  
  # Progress bar
  pb <- txtProgressBar(min = 0,
                       max = length(condi),
                       style = 3)
  
  for (i in seq_len(length(condi))) {
    eCondition <- condi[i]
    # Progress bar
    setTxtProgressBar(pb, i)
    
    # cat("\n\nCONDITION: ",eCondition,"\n##################\n")
    # cat("TECHNICAL REPLICAS\n---------------------------\n")
    # cat("- ", eCondition,"\n")
    
    conditionOne <- data[which(data$Condition == eCondition), ]
    
    # FIRST CHECK FOR TECHNICAL REPLICAS
    bioreplicasAll <- unique(conditionOne$BioReplicate)
    
    for (eBioreplica in bioreplicasAll) {
      # cat('\tChecking for technical replicas in ',eBioreplica, "\n")
      biorepli <-
        conditionOne[conditionOne$BioReplicate == eBioreplica, ]
      here <- unique(biorepli$Run)
      
      if (length(here) > 1) {
        #Limit to 2 technical replicas: (length(here) > 1 & length(here) == 2)
        # We are expecting no more than 2 technical replicas. 
        # If there is more... it is worthy to double check
        # cat('\t\t>>Reproducibility for technical replicas of',eBioreplica,":")
        #Need to change the RUN number to letters (TR: TECHNICAL REPLICA)
        biorepli$TR <- biorepli$Run
        biorepli$TR[biorepli$TR == here[1]] <- 'tr1'
        biorepli$TR[biorepli$TR == here[2]] <- 'tr2'
        
        # Let's select unique features per TECHNICAL REPLICAS, 
        # and sum them up (the same feature might have been many differnt times)
        biorepliaggregated <-
          aggregate(
            Intensity ~ Feature+Proteins+Condition+BioReplicate+Run+TR,
            data = biorepli,
            FUN = sum
          )
        biorepliaggregated$Intensity <-
          log2(biorepliaggregated$Intensity)
        bdc <-
          data.table::dcast(data = biorepliaggregated,
                            Feature + Proteins ~ TR,
                            value.var = 'Intensity')
        
        # Get the number of proteins
        np <- dim(bdc)[1]
        corr_coef <-
          round(cor(bdc$tr1, bdc$tr2, use = "pairwise.complete.obs"),
                digits = 2)
        # cat("r:\t",corr_coef,"\n")
        p1 <- ggplot(bdc, aes(x = tr1, y = tr2))
        p1 <- p1 + geom_point() + geom_rug() + 
          geom_density_2d(colour = 'lightgreen')
        p1 <-
          p1 + geom_smooth(colour = "green",
                           fill = "lightblue",
                           method = 'lm')
        p1 <- p1 + theme_light()
        p1 <-
          p1 + labs(
            title = paste(
              "Reproducibility between Technical Replicas\nBioReplica:",
              eBioreplica,
              "  (n = ",
              np,
              ", r = ",
              corr_coef,
              ")"
            )
          )
        print(p1)
      } else if (length(here) == 1) {
        # cat("\t\tOnly one technical replica\n")
      } else{
        cat(
          "\n>>MORE THAN TWO TECHNICAL REPLICAS IN THIS EXPERIMENTS? That is very strange.\n\n"
        )
        stop("\nPlease, Check the keys files\n")
      }
    } # Checking the reproducibility between Technical Replicas
    
    
    # NOW BETWEEN BIOREPLICAS
    # Before comparing the different biological replicas, 
    # aggregate the technical replicas
    # cat("\nBIOLOGICAL REPLICAS\n---------------------------\n")
    #
    # First choose the maximum for the technical replicas as before, 
    # but first check whether there are more than one
    if (length(here) > 1) {
      conditionOne <-
        aggregate(
          Intensity ~ Feature + Proteins + Condition + BioReplicate + Run,
          data = conditionOne,
          FUN = sum
        )
      b <-
        aggregate(
          Intensity ~ Feature + Proteins + Condition + BioReplicate,
          data = conditionOne,
          FUN = mean
        )
    } else{
      b <-
        aggregate(
          Intensity ~ Feature + Proteins + Condition + BioReplicate,
          data = conditionOne,
          FUN = sum
        )
    }
    
    b$Intensity <- log2(b$Intensity)
    blist <- unique(b$BioReplicate)
    if (length(blist) > 1) {
      # We need at least TWO BIOLOGICAL REPLICAS
      
      to <- length(blist) - 1
      
      for (i in seq_len(to)) {
        j <- i + 1
        for (k in j:length(blist)) {
          br1 <- blist[i]
          br2 <- blist[k]
          
          # cat("\tChecking reproducibility between ",br1, "and ",br2 ,"\t")
          bcfinal <-
            data.table::dcast(data = b,
                              Feature + Proteins ~ BioReplicate,
                              value.var = 'Intensity')
          
          # Let's check the total number of peptides here...
          checkTotalNumber <- subset(bcfinal, select = c(br1, br2))
          # checkTotalNumber <- 
          # checkTotalNumber[complete.cases(checkTotalNumber),]
          
          npt <- dim(checkTotalNumber)[1]
          
          corr_coef <-
            round(cor(bcfinal[[br1]], bcfinal[[br2]], 
                      use = "pairwise.complete.obs"), digits = 2)
          # cat("r:\t",corr_coef,"\n")
          p2 <-
            ggplot(bcfinal, aes(x = bcfinal[[br1]], y = bcfinal[[br2]]))
          p2 <-
            p2 + geom_point()  + geom_rug() + geom_density_2d(colour = 'red')
          p2 <-
            p2 + geom_smooth(colour = "red",
                             fill = "lightgreen",
                             method = 'lm')
          p2 <- p2 + theme_light()
          p2 <-
            p2 + labs(
              title = paste(
                "Peptide Reproducibility between Bioreplicas\n(condition:",
                eCondition,
                ")\n",
                br1,
                "vs",
                br2,
                "(n =",
                npt,
                " r = ",
                corr_coef,
                ")"
              )
            )
          p2 <- p2 + labs(x = br1)
          p2 <- p2 + labs(y = br2)
          print(p2)
        } #end for
      } #end for
    } else{
      cat("\tONLY ONE BIOLOGICAL REPLICA AVAILABLE (plots are not possible)\n")
    }
  } # all the conditions
  # Close Progress bar
  close(pb)
}

# ------------------------------------------------------------------------------
#' @title Plot abundance boxplots
#'
#' @description Plot abundance boxplots
#' @param data (data.frame) processed modelqc
#' @return Abundacen boxplot
#' @keywords internal, plot, abundance
.artms_plotAbundanceBoxplots <- function(data) {
  p1 <-
    ggplot2::ggplot(data, aes(x = SUBJECT_ORIGINAL, 
                              y = ABUNDANCE, 
                              fill = ABUNDANCE))
  p1 <- p1 + geom_boxplot(aes(fill = SUBJECT_ORIGINAL))
  p1 <- p1 + theme_linedraw()
  p1 <-
    p1 + theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
      legend.position = "none"
    )
  p1 <- p1 + labs(x = "BIOREPLICATES")
  p1 <- p1 + ggtitle("Relative Abundance BioReplicates")
  print(p1)
  
  p2 <-
    ggplot2::ggplot(data, aes(
      x = as.factor(GROUP_ORIGINAL),
      y = ABUNDANCE,
      fill = ABUNDANCE
    ))
  p2 <- p2 + geom_boxplot(aes(fill = GROUP_ORIGINAL))
  p2 <- p2 + theme_linedraw()
  p2 <-
    p2 + theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
      legend.position = "none"
    )
  p2 <- p2 + labs(x = "CONDITIONS")
  p2 <- p2 + ggtitle("Relative Abundance Conditions")
  print(p2)
}


# ------------------------------------------------------------------------------
#' @title Total Number of unique proteins based on abundance data
#'
#' @description Total Number of unique proteins per biological replicate and
#' conditions
#' @param data (data.frame) modelqc
#' @return (pdf) Barplots with the number of proteins per br / condition
#' @keywords internal, plots, abundance, counts
.artms_plotNumberProteinsAbundance <- function(data) {
  x <- data[c('PROTEIN', 'SUBJECT_ORIGINAL')]
  y <- unique(x)
  names(y)[grep('SUBJECT_ORIGINAL', names(y))] <- 'BioReplicate'
  z <-
    ggplot2::ggplot(y, aes(x = BioReplicate, fill = BioReplicate))
  z <- z + geom_bar(stat = "count")
  z <-
    z + theme(axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ))
  z <-
    z + geom_text(
      stat = 'count',
      aes(label = ..count..),
      vjust = -0.5,
      size = 2.7
    )
  z <- z + ggtitle("Unique Proteins in BioReplicates")
  print(z)
  
  a <- data[c('PROTEIN', 'GROUP_ORIGINAL')]
  b <- unique(a)
  names(b)[grep('GROUP_ORIGINAL', names(b))] <- 'Condition'
  c <- ggplot2::ggplot(b, aes(x = Condition, fill = Condition))
  c <- c + geom_bar(stat = "count")
  c <-
    c + theme(axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ))
  c <-
    c + geom_text(
      stat = 'count',
      aes(label = ..count..),
      vjust = -0.5,
      size = 2.7
    )
  c <- c + ggtitle("Unique Proteins in Conditions")
  print(c)
}

# ------------------------------------------------------------------------------
#' @title Generate reproducibility plots based on abundance data
#' (normalized intensities from MSstats modelqc)
#'
#' @description Generate reproducibility plots based on abundance data
#' (normalized intensities from MSstats modelqc)
#' @param data (data.frame) Protein abundance (modelqc)
#' @return Reproducibility plots based on abundance data 
#' (normalized intensities)
#' @keywords plot, reproducibility, abundance
.artms_plotReproducibilityAbundance <- function(data) {
  condi <- unique(data$GROUP_ORIGINAL)
  
  # Progress bar
  pb <- txtProgressBar(min = 0,
                       max = length(condi),
                       style = 3)
  
  for (i in seq_len(length(condi))) {
    eCondition <- condi[i]
    
    # Progress bar
    setTxtProgressBar(pb, i)
    
    # cat("\n##################\nCONDITION: ",eCondition,"\n###############\n")
    # cat("TECHNICAL REPLICAS\n---------------------------\n")
    # cat("- ", eCondition,"\n")
    
    conditionOne <- data[which(data$GROUP_ORIGINAL == eCondition), ]
    
    # FIRST CHECK FOR TECHNICAL REPLICAS
    bioreplicasAll <- unique(conditionOne$SUBJECT_ORIGINAL)
    # plot_tr = list()
    for (eBioreplica in bioreplicasAll) {
      # cat('\tChecking for technical replicas in ',eBioreplica, "\n")
      biorepli <-
        conditionOne[conditionOne$SUBJECT_ORIGINAL == eBioreplica, ]
      here <- unique(biorepli$RUN)
      
      if (length(here) > 1) {
        # Check whether we have more than 2 technical replicas and let 
        # the user know:
        if (length(here) > 2) {
          cat(
            "\n\n(-)----- WARNING: More than 2 technical replicas! make sure that this is right\n\n"
          )
        }
        # We are expecting no more than 2 technical replicas. 
        # If there is more... it is worthy to double check
        # cat('\t\t>>Plotting Reproducibility between technical replicas ',eBioreplica,"\n")
        #Need to change the RUN number to letters (TR: TECHNICAL REPLICA)
        biorepli$TR <- biorepli$RUN
        biorepli$TR[biorepli$TR == here[1]] <- 'tr1'
        biorepli$TR[biorepli$TR == here[2]] <- 'tr2'
        bdc <-
          data.table::dcast(data = biorepli, PROTEIN ~ TR, 
                            value.var = 'ABUNDANCE')
        bdc <- bdc[complete.cases(bdc), ]
        # Get the number of proteins
        np <- length(unique(bdc$PROTEIN))
        corr_coef <- round(cor(bdc$tr1, bdc$tr2), digits = 2)
        p1 <- ggplot2::ggplot(bdc, aes(x = tr1, y = tr2))
        p1 <- p1 + geom_point()
        p1 <-
          p1 + geom_smooth(colour = "green",
                           fill = "lightblue",
                           method = 'lm')
        p1 <- p1 + theme_light()
        p1 <-
          p1 + labs(
            title = paste(
              "Reproducibility between Technical Replicas\nBioReplica:",
              eBioreplica,
              "  (n = ",
              np,
              ", r = ",
              corr_coef,
              ")"
            )
          )
        print(p1)
        # plot_tr[[eBioreplica]] <- p1
      } else if (length(here) < 1) {
        cat(
          "\n\n\t(-) ERROR: something is wrong when checking for the number of technical replicas\n\n"
        )
        stop("\nCheck the experiment\n\n")
      }
    } # Checking the reproducibility between Technical Replicas
    
    
    # NOW BETWEEN BIOREPLICAS
    # Before comparing the different biological replicas, 
    # aggregate on the technical replicas
    # cat("\nBIOLOGICAL REPLICAS\n---------------------------\n")
    b <-
      aggregate(
        ABUNDANCE ~ PROTEIN + GROUP_ORIGINAL + SUBJECT_ORIGINAL,
        data = conditionOne,
        FUN = mean
      ) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    blist <- unique(b$SUBJECT_ORIGINAL)
    if (length(blist) > 1) {
      # We need at least TWO BIOLOGICAL REPLICAS
      to <- length(blist) - 1
      for (i in seq_len(to)) {
        j <- i + 1
        for (k in j:length(blist)) {
          br1 <- blist[i]
          br2 <- blist[k]
          # cat("\tChecking reproducibility between ",br1, "and ",br2 ,"\n")
          
          bc <-
            data.table::dcast(data = b,
                              PROTEIN ~ SUBJECT_ORIGINAL,
                              value.var = 'ABUNDANCE')
          bc <- bc[complete.cases(bc), ]
          
          npt <- length(unique(bc$PROTEIN))
          
          corr_coef <- round(cor(bc[[br1]], bc[[br2]]), digits = 2)
          p2 <- ggplot2::ggplot(bc, aes(x = bc[[br1]], y = bc[[br2]]))
          p2 <- p2 + geom_point()
          p2 <-
            p2 + geom_smooth(colour = "red",
                             fill = "lightgreen",
                             method = 'lm')
          p2 <- p2 + theme_light()
          p2 <-
            p2 + labs(
              title = paste(
                "Reproducibility between Bioreplicas (condition:",
                eCondition,
                ")\n",
                br1,
                "vs",
                br2,
                "(n =",
                npt,
                " r = ",
                corr_coef,
                ")"
              )
            )
          p2 <- p2 + labs(x = br1)
          p2 <- p2 + labs(y = br2)
          print(p2)
        }
      }
      cat("\n")
    } else{
      cat("\tONLY ONE BIOLOGICAL REPLICA AVAILABLE (plots are not possible)\n")
    }
  } # all the conditions
  
  close(pb)
}

# ------------------------------------------------------------------------------
#' @title Plot correlation between conditions
#'
#' @description Plot correlation between conditions
#' @param data (data.frame) of Protein Abundance (MSstats modelqc)
#' @param numberBiologicalReplicas (int) Number of biological replicates
#' @return (ggplot.object) A correlation plot between conditions
#' @keywords internal, plot, correlation
.artms_plotCorrelationConditions <-
  function(data, numberBiologicalReplicas) {
    # Before jumping to merging biological replicas:
    # Technical replicas: aggregate on the technical replicas
    b <-
      aggregate(
        ABUNDANCE ~ PROTEIN + GROUP_ORIGINAL + SUBJECT_ORIGINAL,
        data = data,
        FUN = mean
      ) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # Aggregate now the CONDITIONS on the biological replicas:
    
    # One way to do this would be to be very stringent, 
    # requiring to find data in all biological replicas:
    # allBiologicalReplicas <- function(x){
    # ifelse(sum(!is.na(x)) == numberBiologicalReplicas, 
    # mean(x, na.rm = TRUE), NA)}
    # datadc <- data.table::dcast(data=b, 
    # PROTEIN~GROUP_ORIGINAL, value.var = 'ABUNDANCE', 
    # fun.aggregate = allBiologicalReplicas, fill = 0)
    
    # Or a most relaxed way:
    datadc <-
      data.table::dcast(
        data = b,
        PROTEIN ~ GROUP_ORIGINAL,
        value.var = 'ABUNDANCE',
        fun.aggregate = mean
      ) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    # before <- dim(datadc)[1]
    # l <- dim(datadc)[2]
    # datadc <- datadc[apply(datadc[c(2:l)],1,function(z) !any(z==0)),]
    # evenafter <- dim(datadc)[1]
    # datadc <- datadc[complete.cases(datadc),]
    # after <- dim(datadc)[1]
    # cat("Total proteins before: ", before, "\nRemoving the 0s: ",evenafter, "\nTotal proteins (only complete cases): ", after, "\n\n")
    
    blist <- unique(b$GROUP_ORIGINAL)
    if (length(blist) > 1) {
      # We need at least TWO CONDITIONS
      to <- length(blist) - 1
      for (i in seq_len(to)) {
        j <- i + 1
        for (k in j:length(blist)) {
          br1 <- blist[i]
          br2 <- blist[k]
          # cat("\t",br1,"-",br2,":")
          
          npt <- length(unique(datadc$PROTEIN))
          
          corr_coef <-
            round(cor(datadc[[br1]], 
                      datadc[[br2]], 
                      use = "complete.obs"), digits = 2)
          # cat ("r:",corr_coef,"\n")
          
          p2 <-
            ggplot2::ggplot(datadc, aes(x = datadc[[br1]], y = datadc[[br2]]))
          p2 <- p2 + geom_point()
          p2 <-
            p2 + geom_smooth(colour = "blue",
                             fill = "lightblue",
                             method = 'lm')
          p2 <- p2 + theme_light()
          p2 <-
            p2 + labs(
              title = paste(
                "CORRELATION between CONDITIONS:\n",
                br1,
                "and",
                br2,
                "\n(n =",
                npt,
                " r = ",
                corr_coef,
                ")"
              )
            )
          p2 <- p2 + labs(x = br1)
          p2 <- p2 + labs(y = br2)
          print(p2)
        } # FOR loop
      } # For loop
      # cat("\n")
    } else{
      cat("\tONLY ONE BIOLOGICAL REPLICA AVAILABLE (plots are not possible)\n")
    }
  }

# ------------------------------------------------------------------------------
#' @title Plot correlation between quantifications different quantified
#' comparisons
#'
#' @description Plot correlation between all quantifications, i.e., different
#' quantified comparisons
#' @param datai (data.frame) Processed MSstats results
#' @return (ggplot.object) Plot correlation between quantifications and r values
#' @keywords internal, plot, correlation, log2fc
.artms_plotRatioLog2fc <- function(datai) {
  datadc <- dcast(data = datai, Protein ~ Label, value.var = 'log2FC')
  before <- dim(datadc)[1]
  l <- dim(datadc)[2]
  datadc <-
    do.call(data.frame, lapply(datadc, function(x)
      replace(x, is.infinite(x), NA)))
  datadc <- datadc[complete.cases(datadc), ]
  after <- dim(datadc)[1]
  cat(
    "---Total unique identifiers before: ",
    before,
    "\n---Total unique identifiers (only complete cases): ",
    after,
    "\n\n"
  )
  
  blist <- unique(datai$Label)
  blist <- gsub("-", ".", blist)
  
  if (length(blist) > 1) {
    # We need at least TWO CONDITIONS
    to <- length(blist) - 1
    # Progress bar
    pb <- txtProgressBar(min = 0,
                         max = length(blist) - 1,
                         style = 3)
    for (i in seq_len(to)) {
      setTxtProgressBar(pb, i)
      j <- i + 1
      for (k in j:length(blist)) {
        br1 <- blist[i]
        br2 <- blist[k]
        # cat("\tChecking relation between conditions ",br1, "and ",br2 ,":")
        
        npt <- length(unique(datadc$Protein))
        
        corr_coef <-
          round(cor(datadc[[br1]], datadc[[br2]]), digits = 2)
        # cat ("r: ",corr_coef,"\n")
        
        p3 <-
          ggplot2::ggplot(datadc, aes(x = datadc[[br1]], y = datadc[[br2]]))
        p3 <- p3 + geom_point() + geom_rug() + geom_density_2d()
        p3 <-
          p3 + geom_smooth(colour = "red",
                           fill = "lightblue",
                           method = 'lm')
        p3 <- p3 + theme_light()
        p3 <-
          p3 + labs(title = paste0(
            "log2fc(",
            br1,
            ") vs log2fc(",
            br2,
            ")\n(n =",
            npt,
            " r = ",
            corr_coef,
            ")"
          ))
        p3 <- p3 + labs(x = br1)
        p3 <- p3 + labs(y = br2)
        print(p3)
      } # FOR loop
    } # For loop
    close(pb)
  } else{
    cat("\tONLY ONE BIOLOGICAL REPLICA AVAILABLE (plots are not possible)\n")
  }
}

# ------------------------------------------------------------------------------
#' @title Generate PCA plots based on abundance data
#'
#' @description Generate PCA plots based on abundance data
#' @param data Data.frame output from `artms_loadModelQCstrict`
#' @param filename Prefix to generate output names (WITH NO EXTENSION)
#' @param allConditions Conditions selected to generate the plots
#' @return PCA plots based on abundance data (pdf format)
#' @keywords internal, plot, pca
.artms_getPCAplots <- function(data, filename, allConditions) {
  # PRINCIPAL COMPONENT ANALYSIS
  # Using the following packages:
  # FactoMineR, factoextra, corrplot, PerformanceAnalytics
  
  # Remove NA
  nogenename2 <- data[duplicated(data$Gene), ]
  if (dim(nogenename2)[1] > 0) {
    # cat("\t---Removing proteins without a gene name:\n")
    data <- data[complete.cases(data), ]
  }
  
  data <- data[!duplicated(data[, c("Gene")]),]
  rownames(data) <- data$Gene
  df <- data[, allConditions]
  
  # Correlation matrix
  df.cor.matrix <- round(cor(df, use = "complete.obs"), 2)
  
  # Correlation plots:
  out.correlation <- paste0(filename, "-correlations.pdf")
  pdf(out.correlation)
  corrplot::corrplot(
    df.cor.matrix,
    type = "upper",
    tl.col = "black",
    tl.srt = 45,
    diag = FALSE,
    addCoef.col = TRUE
  )
  PerformanceAnalytics::chart.Correlation(df,
                                          histogram = TRUE,
                                          pch = 19,
                                      main = "Correlation between Conditions")
  garbage <- dev.off()
  
  res.pca <-
    FactoMineR::PCA(df,
                    scale.unit = TRUE,
                    ncp = 4,
                    graph = FALSE)
  eigenvalues <- res.pca$eig
  
  out.pca01 <- paste0(filename, "-pca01.pdf")
  pdf(out.pca01)
  par(mfrow = c(1, 1))
  plot(res.pca, choix = "ind", new.plot = FALSE)
  plot(res.pca, choix = "var", new.plot = FALSE)
  # This is equivalent to
  # PCA(df, scale.unit = TRUE, ncp = 4, graph = TRUE)
  garbage <- dev.off()
  
  out.pca02 <- paste0(filename, "-pca02.pdf")
  pdf(out.pca02)
  barplot(
    eigenvalues[, 2],
    names.arg = seq_len(nrow(eigenvalues)),
    main = "Variances",
    xlab = "Principal Components",
    ylab = "Percentage of variances",
    col = "steelblue"
  )
  lines(
    x = seq_len(nrow(eigenvalues)),
    eigenvalues[, 2],
    type = "b",
    pch = 19,
    col = "red"
  )
  garbage <- dev.off()
  
  h <- factoextra::fviz_pca_var(res.pca, 
                                col.var = "contrib") + theme_minimal()
  i <- factoextra::fviz_pca_biplot(res.pca,  
                                   labelsize = 3, 
                                   pointsize = 0.8) + theme_minimal()
  j <- factoextra::fviz_contrib(res.pca, choice = "var", axes = 1)
  l <- factoextra::fviz_contrib(res.pca, choice = "var", axes = 2)
  
  out.pca03 <- paste0(filename, "-pca03.pdf")
  pdf(out.pca03)
  print(h)
  print(i)
  print(j)
  print(l)
  garbage <- dev.off()
}

# ------------------------------------------------------------------------------
#' @title Volcano plot (log2fc / pvalues)
#'
#' @description It generates a scatter-plot used to quickly identify changes
#' @param mss_results (data.frame or file) Selected MSstats results
#' @param lfc_upper (numeric) log2fc upper threshold (positive value)
#' @param lfc_lower (numeric) log2fc lower threshold (negative value)
#' @param whatPvalue (char) `pvalue` or `adj.pvalue` (default)
#' @param FDR (numeric) False Discovery Rate threshold
#' @param output_name (char) Name for the output file
#' @param PDF (logical) Option to generate pdf format. Default: `T`
#' @param decimal_threshold (numeric) Decimal threshold for the pvalue.
#' Default: 16 (10^-16)
#' @keywords plot, volcano
#' @return (pdf) of a volcano plot
#' @examples
#' artms_volcanoPlot(mss_results = artms_data_ph_msstats_results,
#'                   whatPvalue = "pvalue",
#'                   PDF = FALSE)
#' @export
artms_volcanoPlot <- function(mss_results,
                              lfc_upper = 1,
                              lfc_lower = -1,
                              whatPvalue = "adj.pvalue",
                              FDR = 0.05,
                              PDF = TRUE,
                              output_name = '',
                              decimal_threshold = 16) {
  
  cat(">> GENERATING VOLCANO PLOT FROM MSSTATS RESULTS\n")
  if (PDF) {
    if (!grepl("\\.pdf", output_name)) {
      stop("FILE EXTENSION '.pdf' IS MISSED for < output_name >")
    }
  }
  
  mss_results <- .artms_checkIfFile(mss_results)
  
  # handle cases where log2FC is Inf. There are no pvalues or other 
  # information for these cases :(
  # Issues with extreme_val later if we have Inf/-Inf values.
  if (sum(is.infinite(mss_results$log2FC)) > 0) {
    idx <- is.infinite(mss_results$log2FC)
    mss_results$log2FC[idx] <- NA
  }
  
  min_x <- -ceiling(max(abs(mss_results$log2FC), na.rm = TRUE))
  max_x <- ceiling(max(abs(mss_results$log2FC), na.rm = TRUE))
  
  # Deal with special cases in the data where we have pvalues = Inf,NA,0
  if (whatPvalue == "adj.pvalue") {
    if (sum(is.na(mss_results$adj.pvalue)) > 0)
      mss_results <- mss_results[!is.na(mss_results$adj.pvalue), ]
    if (nrow(mss_results[mss_results$adj.pvalue == 0 |
                         mss_results$adj.pvalue == -Inf, ]) > 0)
      mss_results[!is.na(mss_results$adj.pvalue) &
                    (mss_results$adj.pvalue == 0 |
                       mss_results$adj.pvalue == -Inf), ]$adj.pvalue = 10 ^ -decimal_threshold
    max_y = ceiling(-log10(min(mss_results[mss_results$adj.pvalue > 0, ]$adj.pvalue, na.rm = TRUE))) + 1
  } else if (whatPvalue == "pvalue") {
    if (sum(is.na(mss_results$pvalue)) > 0)
      mss_results <- mss_results[!is.na(mss_results$pvalue), ]
    if (nrow(mss_results[mss_results$pvalue == 0 |
                         mss_results$pvalue == -Inf, ]) > 0)
      mss_results[!is.na(mss_results$pvalue) &
                    (mss_results$pvalue == 0 |
                       mss_results$pvalue == -Inf), ]$pvalue = 10 ^ -decimal_threshold
    max_y = ceiling(-log10(min(mss_results[mss_results$pvalue > 0, ]$pvalue, na.rm = TRUE))) + 1
  } else{
    stop("The whatPvalue argument is wrong. Valid options: < pvalue > or < adj.pvalue >")
  }
  l <- length(unique(mss_results$Label))
  w_base <- 7
  h_base <- 7
  
  if (l <= 2) {
    w <- w_base * l
  } else{
    w <- w_base * 2
  }
  h <- h_base * ceiling(l / 2)
  
  if (whatPvalue == "adj.pvalue") {
    p <- ggplot(mss_results, aes(x = log2FC, y = -log10(adj.pvalue)))
    p <- p + geom_point(colour = 'grey') +
      geom_point(
        data = mss_results[mss_results$adj.pvalue <= FDR &
                             mss_results$log2FC >= lfc_upper, ],
        aes(x = log2FC, y = -log10(adj.pvalue)),
        colour = 'red',
        size = 2
      ) +
      geom_point(
        data = mss_results[mss_results$adj.pvalue <= FDR &
                             mss_results$log2FC <= lfc_lower, ],
        aes(x = log2FC, y = -log10(adj.pvalue)),
        colour = 'blue',
        size = 2
      ) +
      geom_vline(xintercept = c(lfc_lower, lfc_upper),
                 lty = 'dashed') +
      geom_hline(yintercept = -log10(FDR), lty = 'dashed') +
      xlim(min_x, max_x) +
      ylim(0, max_y) +
      facet_wrap(facets = ~ Label,
                 ncol = 2,
                 scales = 'fixed')
  } else{
    p <- ggplot(mss_results, aes(x = log2FC, y = -log10(pvalue)))
    p <- p + geom_point(colour = 'grey') +
      geom_point(
        data = mss_results[mss_results$pvalue <= FDR &
                             mss_results$log2FC >= lfc_upper, ],
        aes(x = log2FC, y = -log10(pvalue)),
        colour = 'red',
        size = 2
      ) +
      geom_point(
        data = mss_results[mss_results$pvalue <= FDR &
                             mss_results$log2FC <= lfc_lower, ],
        aes(x = log2FC, y = -log10(pvalue)),
        colour = 'blue',
        size = 2
      ) +
      geom_vline(xintercept = c(lfc_lower, lfc_upper),
                 lty = 'dashed') +
      geom_hline(yintercept = -log10(FDR), lty = 'dashed') +
      xlim(min_x, max_x) +
      ylim(0, max_y) +
      facet_wrap(facets = ~ Label,
                 ncol = 2,
                 scales = 'fixed')
  }
  if (PDF) {
    pdf(output_name, width = w, height = h)
    print(p)
    garbage <- dev.off()
  } else{
    print(p)
  }
}