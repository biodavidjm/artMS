# ------------------------------------------------------------------------------
#' @title Quality Control analysis of the evidence-like metabolomics dataset
#'
#' @description Quality Control analysis of the evidence-like metabolomics 
#' dataset
#' @param evidence_file (char or data.frame) The evidence file path and name, or
#' data.frame
#' @param keys_file (char or data.frame) The keys file path and name or 
#' data.frame
#' @param output_name (char) prefix output name (no extension).
#' Default: "qcPlots_metab"
#' @param met_exp (char) Proteomics experiment. Only one option available
#' (so far):
#' - `MV`: Markview output
#' @param plotINTDIST if `TRUE` (default) plots both *Box-dot plot* 
#' and *Jitter plot* of biological replicates based on MS (raw) 
#' intensity values.
#' @param plotREPRO if `TRUE` (default) plots a correlation dotplot for all the 
#' combinations of biological replicates of conditions, based on MS Intensity 
#' values using features (mz_rt+charge)
#' @param plotCORMAT if `TRUE` (default) generates up to 3 pdf files for 
#' technical replicates, biological replicates, and conditions. Each pdf file 
#' contains: 
#' - *Correlation matrix* for all the biological replicates using 
#' MS Intensity values, 
#' - *Clustering matrix* of the MS Intensities and correlation distribution 
#' - *histogram* of the distribution of correlations
#' @param plotINTMISC if `TRUE` (default) plots several pages, including 
#' bar plots of *Total Sum of Intensities in BioReplicates*, 
#' *Total Sum of Intensities in Conditions*, 
#' *Total Feature Counts in BioReplicates*, 
#' *Total Feature Counts in conditions* separated by categories 
#' (INT: has a intensity value NOINT: no intensity value ) 
#' *Box plots* of MS Intensity values per 
#' biological replicates and conditions; *bar plots* of total intensity 
#' by bioreplicates and conditions; Barplots of 
#' *total feature counts* by bioreplicates and conditions.
#' @param printPDF If `TRUE` (default) prints out the pdfs. Warning: plot
#' objects are not returned due to the large number of them. 
#' @param verbose (logical) `TRUE` (default) shows function messages
#' @return Quality control files and plots for metabolomics
#' @keywords QC, quality, control, evidence, metabolomics
#' @examples
#' # Testing that input arguments cannot be null
#' artmsQualityControlMetabolomics(evidence_file = NULL,
#'                  keys_file = NULL,
#'                  met_exp = "MV")
#' @export
artmsQualityControlMetabolomics <- function(evidence_file,
                             keys_file,
                             met_exp = c('MV'),
                             output_name = "qcPlots_metab",
                             plotINTDIST = TRUE,
                             plotREPRO = TRUE,
                             plotCORMAT = TRUE,
                             plotINTMISC = TRUE,
                             printPDF = TRUE,
                             verbose = TRUE) {
  
  if(any(missing(evidence_file) | missing(keys_file)))
    stop("Missed (one or many) required argument(s)
         Please, check the help of this function to find out more")
  
  if( is.null(evidence_file) & is.null(keys_file) ){
    return("Both <evidence_file> and <keys_file> cannot be NULL")
  }

  met_exp <- toupper(met_exp)
  met_exp <- match.arg(met_exp)
  
  if(verbose) message(" QUALITY CONTROL -------------------\n>> LOADING FILES ")
  
  # EVIDENCE:
  evidencekeys <- artmsMergeEvidenceAndKeys(evidence_file, 
                                             keys_file,
                                             verbose = verbose)
  
  # remove 
  
  ekselecta <- evidencekeys
  ekselectaBioreplica <- evidencekeys
  
  # Checking the overall distribution of intensities before anything else
  # Based on Intensity
  ekselectaBioreplica$Abundance <- log2(ekselectaBioreplica$Intensity)
  
  if (plotINTDIST){
    if(verbose) message(">> GENERATING THE INTENSITY DISTRIBUTION PLOTS ")
    intDistribution <- paste0(output_name, ".qcplot.plotINTDIST.pdf")
  
    j <- ggplot2::ggplot(ekselectaBioreplica, aes(BioReplicate, Intensity))
    j <- j + geom_jitter(width = 0.3, size = 0.5)
    j <- j + theme_minimal()
    j <-
      j + theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ))
    # Based on Abundance
    k <- ggplot2::ggplot(ekselectaBioreplica, aes(BioReplicate, Abundance))
    k <- k + geom_jitter(width = 0.3, size = 0.5)
    k <- k + theme_minimal()
    k <-
      k + theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ))

    if(printPDF) pdf(intDistribution)
    plot(j)
    plot(k)
    if(printPDF) garbage <- dev.off()
    
  }
  
  # Feature generation: Combine Sequence and Charge.
  evidencekeys$Feature <- 
    paste0(evidencekeys$Modified.sequence, "_", evidencekeys$Charge)
  
  # ============================================================================
  # GENERAL QUALITY CONTROL: CHECK PROPORTION OF IntDetectionS
  # sum of total intensities, label them as IntDetections and non-IntDetections
  # and plot intensities for each group

  evigeneral <-
    evidencekeys[c(
      'Feature',
      'Proteins',
      'Intensity',
      'Condition',
      'BioReplicate',
      'Run'
    )]
  
  # Now let's classify the proteins as INT and INT
  evigeneral$IntDetection <- ifelse( evigeneral$Intensity > 0, "INT", "NOINT")
  
  # AGGREGATE ALL THE INTENSITIES PER PROTEIN, summing everything up
  evigeneral$Intensity <- as.numeric(evigeneral$Intensity)
  
  evisummary <- aggregate(
      Intensity ~ Condition + BioReplicate + IntDetection + Run,
      data = evigeneral,
      FUN = sum
    )
  
  # CLEANING THE EVIDENCE OF NOINT
  if( any( dim(evidencekeys[which(evidencekeys$Intensity == 0),])[1] > 0 ) )
    evidencekeysclean <- evidencekeys[-which(evidencekeys$Intensity == 0),]
  
  if(plotREPRO){
    if(verbose) message(">> GENERATING THE REPRODUCIBILITY PLOTS 
      (Warning: it might take some time) ")
    seqReproName <-
      paste0(output_name, ".qcplot.plotREPRO.pdf")
    
    if(printPDF) pdf(seqReproName)
    .artms_plotReproducibilityEvidence(evidencekeysclean, verbose = verbose)
    if(printPDF) garbage <- dev.off()
  }
  
  # Create matrix of reproducibility TECHNICAL REPLICAS
  data2matrix <- evidencekeysclean
  # Make sure that the Intensity is numeric
  data2matrix$Intensity <- as.numeric(data2matrix$Intensity)
  
  # Check the number of TECHNICAL REPLICAS by 
  # checking the first technical replica
  technicalReplicas <-
    unique(data2matrix$Run[which(data2matrix$BioReplicate == data2matrix$BioReplicate[1])])
  palette.breaks <- seq(1, 3, 0.1)
  color.palette  <-
    colorRampPalette(c("white", "steelblue"))(length(palette.breaks))
  
  if(plotCORMAT){
    if(verbose) message(">> GENERATING CORRELATION MATRICES ")
    if (length(technicalReplicas) > 1) {
      # First aggregate at the protein level by summing up everything
      biorepliaggregated <-
        aggregate(
          Intensity ~ Feature + Proteins + Condition + BioReplicate + Run,
          data = data2matrix,
          FUN = sum
        )
      biorepliaggregated$Intensity <-
        log2(biorepliaggregated$Intensity)
      evidencekeyscleanDCASTbioreplicas <-
        data.table::dcast(data = biorepliaggregated,
                          Proteins + Feature ~ BioReplicate + Run,
                          value.var = "Intensity")
      precordfBioreplicas <-
        evidencekeyscleanDCASTbioreplicas[, 3:dim(evidencekeyscleanDCASTbioreplicas)[2]]
      Mtechnicalrep <-
        cor(precordfBioreplicas, use = "pairwise.complete.obs")
      
      theTechCorDis <-
        .artms_plotCorrelationDistribution(Mtechnicalrep)
      
      # And now for clustering
      if(verbose) message("--- By Technical replicates ")
      matrixCorrelationBioreplicas <-
        paste0(output_name, ".qcplot.plotCORMAT.pdf")
      
      if(printPDF) pdf(matrixCorrelationBioreplicas, width = 20, height = 20) #
      corrplot::corrplot(
        Mtechnicalrep,
        method = "square",
        addCoef.col = "white",
        number.cex = 1,
        tl.cex = 1.5,
        tl.col = "black",
        title = "Matrix Correlation based on Feature Intensities"
      )
      pheatmap::pheatmap(
        Mtechnicalrep,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        cellheight = 10,
        cellwidth = 25,
        main = "Clustering Technical Replicates",
        fontsize = 6,
        fontsize_row = 8,
        fontsize_col = 12,
        border_color = 'black',
        fontfamily = "Helvetica",
        treeheight_row = FALSE,
        treeheight_col = FALSE,
        color = color.palette
      )
      print(theTechCorDis)
      if(printPDF) garbage <- dev.off()
    } else{
      if(verbose) message("--- NO Technical Replicates detected ")
    }
    
    # biological replicates
    biorepliaggregated <- NULL
    
    # First deal with the technical replica
    if (length(technicalReplicas) > 1) {
      # Aggregate at the protein level by summing up everything
      biorepliaggregated <-
        aggregate(
          Intensity ~ Feature + Proteins + Condition + BioReplicate + Run,
          data = data2matrix,
          FUN = sum
        )
      biorepliaggregated <-
        aggregate(
          Intensity ~ Feature + Proteins + Condition + BioReplicate,
          data = biorepliaggregated,
          FUN = mean
        )
    } else{
      biorepliaggregated <-
        aggregate(
          Intensity ~ Feature + Proteins + Condition + BioReplicate,
          data = data2matrix,
          FUN = sum
        )
    }
    
    biorepliaggregated$Intensity <- log2(biorepliaggregated$Intensity)
    evidencekeyscleanDCASTbioreplicas <-
      data.table::dcast(data = biorepliaggregated,
                        Proteins + Feature ~ BioReplicate,
                        value.var = "Intensity")
    precordfBioreplicas <-
      evidencekeyscleanDCASTbioreplicas[, 3:dim(evidencekeyscleanDCASTbioreplicas)[2]]
    Mbioreplicas <-
      cor(precordfBioreplicas, use = "pairwise.complete.obs")
    
    theBiorCorDis <- .artms_plotCorrelationDistribution(Mbioreplicas)
    
    if(verbose) message("--- By Biological replicates ")
    matrixCorrelationBioreplicas <-
      paste0(output_name, ".qcplot.correlationMatrixBR.pdf")
    if(printPDF) pdf(matrixCorrelationBioreplicas,
        width = 20,
        height = 20) #, width = 20, height = 20
    corrplot::corrplot(
      Mbioreplicas,
      method = "square",
      addCoef.col = "white",
      number.cex = 0.9,
      tl.cex = 1.5,
      tl.col = "black",
      title = "Matrix Correlation based on Feature Intensities"
    )
    pheatmap::pheatmap(
      Mbioreplicas,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      cellheight = 10,
      cellwidth = 25,
      main = "Clustering Biological Replicates",
      fontsize = 6,
      fontsize_row = 8,
      fontsize_col = 12,
      border_color = 'black',
      fontfamily = "Helvetica",
      treeheight_row = FALSE,
      treeheight_col = FALSE,
      color = color.palette
    )
    print(theBiorCorDis)
    if(printPDF) garbage <- dev.off()
    
    # Create matrix of reproducibility CONDITIONS
    biorepliaggregated <- NULL
    
    # biological replicas
    # First deal with the technical replica
    if (length(technicalReplicas) > 1) {
      biorepliaggregated <-
        aggregate(
          Intensity ~ Feature + Proteins + Condition + BioReplicate + Run,
          data = data2matrix,
          FUN = sum
        )
      biorepliaggregated <-
        aggregate(
          Intensity ~ Feature + Proteins + Condition + BioReplicate,
          data = biorepliaggregated,
          FUN = mean
        )
    } else{
      biorepliaggregated <-
        aggregate(
          Intensity ~ Feature + Proteins + Condition + BioReplicate,
          data = data2matrix,
          FUN = sum
        )
    }
    # Now let's sum up based on conditions
    biorepliaggregated <-
      aggregate(Intensity ~ Feature + Proteins + Condition,
                data = biorepliaggregated,
                FUN = median)
    biorepliaggregated$Intensity <- log2(biorepliaggregated$Intensity)
    evidencekeyscleanDCASTconditions <-
      data.table::dcast(data = biorepliaggregated,
                        Proteins + Feature ~ Condition,
                        value.var = "Intensity")
    precordfConditions <-
      evidencekeyscleanDCASTconditions[, 3:dim(evidencekeyscleanDCASTconditions)[2]]
    Mcond <- cor(precordfConditions, use = "pairwise.complete.obs")
    
    theCondCorDis <- .artms_plotCorrelationDistribution(Mcond)
    
    if(verbose) message("--- By Conditions ")
    matrixCorrelationCond <-
      paste0(output_name, ".qcplot.correlationMatrixConditions.pdf")
    if(printPDF) pdf(matrixCorrelationCond)
    corrplot(
      Mcond,
      method = "square",
      addCoef.col = "white",
      number.cex = 0.6,
      tl.col = "black",
      title = "Matrix Correlation based on Feature Intensities"
    )
    pheatmap(
      Mcond,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      cellheight = 10,
      cellwidth = 25,
      main = "Clustering Conditions",
      fontsize = 6,
      fontsize_row = 8,
      fontsize_col = 12,
      border_color = 'black',
      fontfamily = "Helvetica",
      treeheight_row = FALSE,
      treeheight_col = FALSE,
      color = color.palette
    )
    print(theCondCorDis)
    if(printPDF) garbage <- dev.off()
  } # END OF plotCORMAT

  
  if(plotINTMISC){
    # DETAILS
    if(verbose) message(">> GENERATING INTENSITY STATS PLOTS ")
    if (met_exp == "MV") {
      ekselect <- evidencekeysclean[c('Feature',
                                      'Proteins',
                                      'Intensity',
                                      'Condition',
                                      'BioReplicate',
                                      'Run')]
      # Aggregate the technical replicas
      ekselecta <-
        aggregate(
          Intensity ~ Proteins + Condition + BioReplicate + Run,
          data = ekselect,
          FUN = sum
        )
      ekselectaBioreplica <-
        aggregate(
          Intensity ~ Proteins + Condition + BioReplicate + Run,
          data = ekselecta,
          FUN = sum
        )
      # Select Unique number of proteins per Condition
      ac <- ekselecta[c('Proteins', 'Condition')]
      b <- unique(ac)
      # Select Unique number of proteins per Biological Replica
      bc <- ekselecta[c('Proteins', 'Condition', 'BioReplicate')]
      bb <- unique(bc)
      # Select unique number of proteins in Technical Replicates
      cc <- ekselecta[c('Proteins', 'Condition', 'BioReplicate', 'Run')]
      cc$TR <- paste0(ekselecta$BioReplicate, "_", ekselecta$Run)
      ccc <- cc[c('Proteins', 'Condition', 'TR')]
      cd <- unique(ccc)
    } else {
      stop("METABOLOMICS experiment not recognized ")
    }
    
    if(verbose) message("--- ", met_exp, " PROCESSED ")
    reproName <- paste0(output_name, ".qcplot.plotINTMISC.pdf")
    
    #QC: SUM of intensities per biological replica (peptides vs IntDetection)
    pisa <-
      ggplot2::ggplot(evisummary,
                      aes(x = BioReplicate, y = Intensity, fill = IntDetection)) +
      geom_bar(stat = "identity",
               position = position_dodge(width = 0.7),
               width = 0.7) +
      theme_minimal() +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
      legend.title = element_blank()) +
      ggtitle("QC: Total Sum of Intensities in BioReplicates") +
      scale_fill_brewer(palette = "Paired")
    
    #QC: SUM of intensities per condition (peptides vs IntDetection)
    pisb <-
      ggplot2::ggplot(evisummary, aes(x = Condition, y = Intensity, 
                                      fill = IntDetection)) +
      geom_bar(stat = "identity",
               position = position_dodge(width = 0.7),
               width = 0.7) +
      theme_minimal() +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
      legend.title = element_blank()) +
      ggtitle("QC: Total Sum of Intensities in Conditions") +
      scale_fill_brewer(palette = "Paired")
    
    #QC: TOTAL COUNT OF PEPTIDES IN EACH BIOLOGICAL REPLICA
    pisc <-
      ggplot2::ggplot(evigeneral, aes(x = BioReplicate, fill = IntDetection)) +
      geom_bar(stat = "count",
               position = position_dodge(width = 0.7),
               width = 0.7) +
      theme_minimal() +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
      legend.title = element_blank()) +
      ggtitle("QC: Total Feature Counts in BioReplicates")
    
    #QC: TOTAL COUNT OF PEPTIDES IN EACH BIOLOGICAL REPLICA
    pisd <-
      ggplot2::ggplot(evigeneral, aes(x = Condition, fill = IntDetection)) +
      geom_bar(stat = "count",
               position = position_dodge(width = 0.7),
               width = 0.7) +
      theme_minimal() +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
      legend.title = element_blank()) +
      ggtitle("QC: Feature Counts in Conditions")
    
    #QC: MAX INTENSITY IN REDUNDANT PEPTIDES, AND AGGREGATE FOR EACH PROTEIN 
    #THE SUM OF INTENSITY
    ekselectaBioreplica2plot <- ekselectaBioreplica
    ekselectaBioreplica2plot$Intensity <-
      log2(ekselectaBioreplica2plot$Intensity)
    pise <-
      ggplot2::ggplot(ekselectaBioreplica2plot,
                      aes(
                        x = as.factor(BioReplicate),
                        y = Intensity,
                        fill = Intensity
                      )) +
      geom_boxplot(aes(fill = BioReplicate)) +
      # scale_y_log10() + coord_trans(y = "log10") +
      theme_linedraw() +
      theme(
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5
        ),
        legend.position = "none"
      ) +
      labs(x = "BioReplicate", y = "log2(Intensity)") +
      ggtitle("Feature Intensity in BioReplicates Excluding IntDetections. 
              Max intensity of TR")
    
    pisf <-
      ggplot2::ggplot(ekselectaBioreplica2plot,
                      aes(x = Condition, y = Intensity, fill = Intensity)) +
      geom_boxplot(aes(fill = Condition)) +
      # scale_y_log10() + coord_trans(y = "log10") +
      theme_linedraw() +
      theme(
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5
        ),
        legend.position = "none"
      ) +
      labs(x = "Condition", y = "log2(Intensity)") +
      ggtitle("Feature Intensity in Conditions Excluding IntDetections. 
              Max intensity of TR")
    
    pisg <- ggplot2::ggplot(ekselectaBioreplica) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      )) +
      geom_bar(
        aes(BioReplicate, Intensity),
        position = "dodge",
        stat = "summary",
        fun.y = "mean",
        fill = "black",
        colour = "orange"
      ) +
      ggtitle("Total Intensity in Biological Replicas Excluding IntDetections. 
              Max intensity of TR")
    
    pish <- ggplot2::ggplot(ekselectaBioreplica) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      )) +
      geom_bar(
        aes(Condition, Intensity),
        position = "dodge",
        stat = "summary",
        fun.y = "mean",
        fill = "black",
        colour = "green"
      ) +
      ggtitle("Total Intensity in Conditions Excluding IntDetections. 
              Max intensity of TR")
    
    pisi <- ggplot2::ggplot(cd, aes(x = TR, fill = Condition)) +
      geom_bar(stat = "count") +
      theme(
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.9
        ),
        legend.position = "none"
      ) +
      geom_text(
        stat = 'count',
        aes(label = ..count..),
        vjust = -0.5,
        size = 2.7
      ) +
      ggtitle("Unique Features in Technical Replicas")
    
    pisj <- ggplot2::ggplot(bb, aes(x = BioReplicate, fill = Condition)) +
      geom_bar(stat = "count") +
      theme(
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5
        ),
        legend.position = "none"
      ) +
      geom_text(
        stat = 'count',
        aes(label = ..count..),
        vjust = -0.5,
        size = 2.7
      ) +
      ggtitle("Unique Features in Biological Replicas")
    
    pisk <- ggplot2::ggplot(b, aes(x = Condition, fill = Condition)) +
      geom_bar(stat = "count") +
      theme(
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5
        ),
        legend.position = "none"
      ) +
      geom_text(
        stat = 'count',
        aes(label = ..count..),
        vjust = -0.5,
        size = 2.7
      ) +
      ggtitle("Unique Features in Condition")
    
    if(printPDF) pdf(reproName)
    print(pisa)
    print(pisb)
    print(pisc)
    print(pisd)
    print(pise)
    print(pisf)
    print(pisg)
    print(pish)
    print(pisi)
    print(pisj)
    print(pisk)
    if(printPDF) garbage <- dev.off()
  } #END OF plotINTMISC
  
  if(verbose) message(">> BASIC QUALITY CONTROL ANALYSIS COMPLETED!  ")
}
