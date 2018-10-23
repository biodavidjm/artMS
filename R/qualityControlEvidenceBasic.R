# ------------------------------------------------------------------------------
#' @title Quality Control analysis of the MaxQuant evidence file
#'
#' @description Quality Control analysis of the MaxQuant evidence file
#' @param evidence_file (char or data.frame) The evidence file path and name, or
#' data.frame
#' @param keys_file (char or data.frame) The keys file path and name or 
#' data.frame
#' @param output_name (char) prefix output name (no extension).
#' Default: "qcPlots_evidence"
#' @param prot_exp (char) Proteomics experiment. 4 options available:
#' - `APMS`: affinity purification mass spectrometry
#' - `AB`: protein abundance
#' - `PH`: protein phosphorylation
#' - `UB`: protein ubiquitination (aka ubiquitylation)
#' @param fractions (binary) Is a fractionated experiment?
#' - 1 yes
#' - 0 no (default)
#' @param plotINTDIST if `TRUE` (default) plots both *Box-dot plot* 
#' and *Jitter plot* of biological replicates based on MS (raw) 
#' intensity values.
#' @param plotREPRO if `TRUE` (default) plots a correlation dotplot for all the 
#' combinations of biological replicates of conditions, based on MS Intensity 
#' values using features (peptide+charge)
#' @param plotCORMAT if `TRUE` (default) plots a 
#' - *Correlation matrix* for all the biological replicates using 
#' MS Intensity values, 
#' - *Clustering matrix* of the MS Intensities and correlation distribution 
#' - *histogram* of the distribution of correlations
#' @param plotINTMISC if `TRUE` (default) plots several pages, including 
#' bar plots of *Total Sum of Intensities in BioReplicates*, 
#' *Total Sum of Intensities in Conditions*, 
#' *Total Peptide Counts in BioReplicates*, 
#' *Total Peptide Counts in conditions* separated by categories: 
#' `CON`: contaminants, `PROT` peptides, `REV` reversed sequences used by 
#' MaxQuant to estimate the FDR; *Box plots* of MS Intensity values per 
#' biological replicates and conditions; *bar plots* of total intensity 
#' (excluding contaminants) by bioreplicates and conditions; Barplots of 
#' *total feature counts* by bioreplicates and conditions.
#' @param plotPTMSTATS IF `TRUE` (default) plots stats related to the 
#' selected modification, including: 
#' *bar plot of peptide counts and intensities*, broken by `PTM/other` 
#' categories; bar plots of *total sum-up of MS intensity values* by 
#' other/PTM categories.
#' @param printPDF If `TRUE` (default) prints out the pdfs. Warning: plot
#' objects are not returned due to the large number of them. 
#' @param verbose (logical) `TRUE` (default) shows function messages
#' @return Quality control files and plots
#' @keywords QC, quality, control, evidence
#' @examples
#' artmsQualityControlEvidenceBasic(evidence_file = artms_data_ph_evidence,
#'                                  keys_file = artms_data_ph_keys, 
#'                                  prot_exp =  "PH", 
#'                                  plotINTDIST = FALSE,
#'                                  plotREPRO = TRUE,
#'                                  plotCORMAT = FALSE,
#'                                  plotINTMISC = FALSE,
#'                                  plotPTMSTATS = FALSE,
#'                                  printPDF = FALSE,
#'                                  verbose = FALSE)
#' 
#' # But we recommend the following test:
#' # 1. Go to a working directory: 
#' # setwd("/path/to/your/working/directory/")
#' # 2. Run the following command to print out all the pdf files
#' # artmsQualityControlEvidenceBasic(evidence_file = artms_data_ph_evidence,
#' #                                  keys_file = artms_data_ph_keys, 
#' #                                  prot_exp =  "PH")
#' # 3. Check your working directory and you should find pdf files with 
#' # all the QC plots
#' @export
artmsQualityControlEvidenceBasic <- function(evidence_file,
                             keys_file,
                             prot_exp = c('AB', 'PH', 'UB', 'APMS'),
                             fractions = 0,
                             output_name = "qcPlots_evidence",
                             plotINTDIST = TRUE,
                             plotREPRO = TRUE,
                             plotCORMAT = TRUE,
                             plotINTMISC = TRUE,
                             plotPTMSTATS = TRUE,
                             printPDF = TRUE,
                             verbose = TRUE) {
  
  if(any(missing(evidence_file) | missing(keys_file)))
    stop("Missed (one or many) required argument(s)
         Please, check the help of this function to find out more")

  prot_exp <- toupper(prot_exp)
  prot_exp <- match.arg(prot_exp)
  supportedExperiments <- c('AB', 'PH', 'UB', 'APMS')
  
  if (any(!prot_exp %in% supportedExperiments)) {
    stop(prot_exp, " is currently not supported.
The experiments supported are: ",
         sprintf('\t%s\n', supportedExperiments))
  }
  
  if (fractions) {
    # Check that the keys file is correct
    keys <- .artms_checkIfFile(keys_file)
    keys <- .artms_checkRawFileColumnName(keys)
    if (any(!'FractionKey' %in% colnames(keys))) {
      stop(' <fractions> was activated but <fractionkey> column not found in the keys file ')
    }
  }
  
  if(verbose) message(" QUALITY CONTROL -------------------\n>> LOADING FILES ")
  
  # EVIDENCE:
  evidencekeys <- artmsMergeEvidenceAndKeys(evidence_file, 
                                             keys_file,
                                             verbose = verbose)
  
  ekselecta <-
    aggregate(Intensity ~ Proteins + Condition + BioReplicate + Run,
              data = evidencekeys,
              FUN = sum)
  ekselectaBioreplica <-
    aggregate(Intensity ~ Proteins + Condition + BioReplicate,
              data = ekselecta,
              FUN = sum)
  
  # Checking the overall distribution of intensities before anything else
  # Based on Intensity
  ekselectaBioreplica$Abundance <-
    log2(ekselectaBioreplica$Intensity)
  
  if(plotINTDIST){
    if(verbose) message(">> GENERATING THE INTENSITY DISTRIBUTION PLOTS ")
    intDistribution <-
      paste0(output_name, ".qcplot.IntensityDistributions.pdf")
    
    j <- ggplot(ekselectaBioreplica, aes(BioReplicate, Intensity))
    j <- j + geom_jitter(width = 0.3, size = 0.5)
    j <- j + theme_minimal()
    j <-
      j + theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ))
    # Based on Abundance
    k <- ggplot(ekselectaBioreplica, aes(BioReplicate, Abundance))
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
  # GENERAL QUALITY CONTROL: CHECK PROPORTION OF CONTAMINANTS
  # sum of total intensities, label them as contaminants and non-contaminants
  # and plot intensities for each group
  
  # Careful with old versions of MaxQuant:
  if (any(grep("Leading.Proteins", names(evidencekeys)))) {
    evidencekeys <-
      artmsChangeColumnName(evidencekeys, 
                             "Leading.Proteins", 
                             "Leading.proteins")
  }
  
  # Combine all the fractions if this is a fractioning experiment by summing
  # them up
  if (fractions) {
    # Sum up all the fractions first
    evidencekeys <-
      aggregate(
        Intensity~Feature+Proteins+Leading.proteins+Condition+BioReplicate +
          Run,
        data = evidencekeys,
        FUN = sum
      )
  }
  
  evigeneral <-
    evidencekeys[c(
      'Feature',
      'Proteins',
      'Leading.proteins',
      'Intensity',
      'Condition',
      'BioReplicate',
      'Run'
    )]
  
  # Now let's classify the proteins as contaminants and no contaminants
  evigeneral$Contaminant <-
    ifelse(
      grepl("CON_|REV_", evigeneral$Leading.proteins),
      ifelse(grepl("REV_", evigeneral$Leading.proteins), "REV", "CON"),
      "PROT"
    )
  
  # AGGREGATE ALL THE INTENSITIES PER PROTEIN, summing everything up
  evigeneral$Intensity <- as.numeric(evigeneral$Intensity)
  
  evisummary <-
    aggregate(
      Intensity ~ Condition + BioReplicate + Contaminant + Run,
      data = evigeneral,
      FUN = sum
    )
  
  # CLEANING THE EVIDENCE OF CONTAMINANTS
  evidencekeysclean <-
    artmsFilterEvidenceContaminants(x = evidencekeys, verbose = verbose)
  
  if (prot_exp == "UB") {
    evidencekeysclean <-
      evidencekeysclean[c(
        'Feature',
        'Modified.sequence',
        'Proteins',
        'Intensity',
        'Condition',
        'BioReplicate',
        'Run'
      )]
    evidencekeysclean <-
      evidencekeysclean[grep("(gl)", evidencekeysclean$Modified.sequence), ]
  } else if (prot_exp == "PH") {
    evidencekeysclean <-
      evidencekeysclean[c(
        'Feature',
        'Modified.sequence',
        'Proteins',
        'Intensity',
        'Condition',
        'BioReplicate',
        'Run'
      )]
    evidencekeysclean <-
      evidencekeysclean[grep("(ph)", evidencekeysclean$Modified.sequence), ]
  }
  
  if(plotREPRO){
    if(verbose) message(">> GENERATING THE REPRODUCIBILITY PLOTS 
      (Warning: it might take some time) ")
    seqReproName <-
      paste0(output_name, ".qcplot.basicReproducibility.pdf")
    
    
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
    if(verbose) message(">> GENERATING CORRELATION MATRICES")
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
        paste0(output_name, ".qcplot.correlationMatrixTR.pdf")
      
      if(printPDF) pdf(matrixCorrelationBioreplicas,
          width = 20,
          height = 20) #
      corrplot(
        Mtechnicalrep,
        method = "square",
        addCoef.col = "white",
        number.cex = 1,
        tl.cex = 1.5,
        tl.col = "black",
        title = "Matrix Correlation based on peptide Intensities"
      )
      pheatmap(
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
    if(printPDF) pdf(matrixCorrelationBioreplicas, width = 20, height = 20)
      corrplot(
        Mbioreplicas,
        method = "square",
        addCoef.col = "white",
        number.cex = 0.9,
        tl.cex = 1.5,
        tl.col = "black",
        title = "Matrix Correlation based on Protein Intensities"
      )
      pheatmap(
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
        title = "Matrix Correlation based on protein Intensities"
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
  } # plotCORMAT ends

  
  if(plotINTMISC){
    # DETAILS
    if(verbose) message(">> GENERATING INTENSITY STATS PLOTS ")
    if (prot_exp == "APMS" | prot_exp == "AB") {
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
    } else if (prot_exp == "UB") {
      ekselectall <-
        evidencekeysclean[c(
          'Feature',
          'Modified.sequence',
          'Proteins',
          'Intensity',
          'Condition',
          'BioReplicate',
          'Run'
        )]
      ekselectgly <-
        ekselectall[grep("(gl)", ekselectall$Modified.sequence), ]
      # Select the REDUNDANT peptides with the meanimum intensity
      ekselectaBioreplica <-
        aggregate(
          Intensity ~ Proteins + Condition + BioReplicate + Run,
          data = ekselectgly,
          FUN = sum
        )
      # Select Unique number of proteins per Condition
      ac <- ekselectaBioreplica[c('Proteins', 'Condition')]
      b <- unique(ac)
      # Select Unique number of proteins per Biological Replica
      bc <-
        ekselectaBioreplica[c('Proteins', 'Condition', 'BioReplicate')]
      bb <- unique(bc)
      # Select unique number of proteins in Technical Replicates
      cc <-
        ekselectaBioreplica[c('Proteins', 'Condition', 'BioReplicate', 'Run')]
      cc$TR <-
        paste0(ekselectaBioreplica$BioReplicate,
               "_",
               ekselectaBioreplica$Run)
      ccc <- cc[c('Proteins', 'Condition', 'TR')]
      cd <- unique(ccc)
      
      # Check the total number of modified and non-modified residues to plot
      evidencekeys$MODIFICATION <-
        ifelse(grepl("(gl)", evidencekeys$Modified.sequence),
               "ub",
               "other")
    } else if (prot_exp == "PH") {
      ekselectall <-
        evidencekeysclean[c(
          'Feature',
          'Modified.sequence',
          'Proteins',
          'Intensity',
          'Condition',
          'BioReplicate',
          'Run'
        )]
      ekselectgly <-
        ekselectall[grep("(ph)", ekselectall$Modified.sequence), ]
      # SUM UP all the peptide intensities for all the proteins: 
      # only one intensity value per protein and biological data
      ekselectaBioreplica <-
        aggregate(
          Intensity ~ Proteins + Condition + BioReplicate + Run,
          data = ekselectgly,
          FUN = sum
        )
      # Select Unique number of proteins per Condition
      ac <- ekselectaBioreplica[c('Proteins', 'Condition')]
      b <- unique(ac)
      # Select Unique number of proteins per Biological Replica
      bc <-
        ekselectaBioreplica[c('Proteins', 'Condition', 'BioReplicate')]
      bb <- unique(bc)
      # Select unique number of proteins in Technical Replicates
      cc <-
        ekselectaBioreplica[c('Proteins', 'Condition', 'BioReplicate', 'Run')]
      cc$TR <-
        paste0(ekselectaBioreplica$BioReplicate,
               "_",
               ekselectaBioreplica$Run)
      ccc <- cc[c('Proteins', 'Condition', 'TR')]
      cd <- unique(ccc)
      
      # Check the total number of modified and non-modified residues to plot
      evidencekeys$MODIFICATION <-
        ifelse(grepl("(ph)", evidencekeys$Modified.sequence),
               "ph",
               "other")
    } else {
      stop("Proteomics experiment not recognized ")
    }
    
    #QC: SUM of intensities per biological replica (peptides vs contaminant)
    pisa <-
      ggplot(evisummary,
             aes(x = BioReplicate, y = Intensity, fill = Contaminant)) +
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
    
    #QC: SUM of intensities per condition (peptides vs contaminant)
    pisb <-
      ggplot(evisummary, aes(x = Condition, y = Intensity, fill = Contaminant)) +
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
      ggplot(evigeneral, aes(x = BioReplicate, fill = Contaminant)) +
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
      ggtitle("QC: Total Peptide Counts in BioReplicates")
    
    #QC: TOTAL COUNT OF PEPTIDES IN EACH BIOLOGICAL REPLICA
    pisd <-
      ggplot(evigeneral, aes(x = Condition, fill = Contaminant)) +
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
      ggtitle("QC: Peptide Counts in Conditions")
    
    #QC: MAX INTENSITY IN REDUNDANT PEPTIDES, AND AGGREGATE FOR EACH PROTEIN 
    #THE SUM OF INTENSITY
    ekselectaBioreplica2plot <- ekselectaBioreplica
    ekselectaBioreplica2plot$Intensity <-
      log2(ekselectaBioreplica2plot$Intensity)
    pise <-
      ggplot(ekselectaBioreplica2plot,
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
      ggtitle("Protein Intensity in BioReplicates\nExcluding contaminants. 
              Max intensity of TR")
    
    pisf <-
      ggplot(ekselectaBioreplica2plot,
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
      ggtitle("Protein Intensity in Conditions\nExcluding contaminants. 
              Max intensity of TR")
    
    pisg <- ggplot(ekselectaBioreplica) +
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
      ggtitle("Total Intensity in Biological Replicas\nExcluding contaminants. 
              Max intensity of TR")
    
    pish <- ggplot(ekselectaBioreplica) +
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
      ggtitle("Total Intensity in Conditions Excluding contaminants. 
              Max intensity of TR")
    
    pisi <- ggplot(cd, aes(x = TR, fill = Condition)) +
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
    
    pisj <- ggplot(bb, aes(x = BioReplicate, fill = Condition)) +
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
    
    pisk <- ggplot(b, aes(x = Condition, fill = Condition)) +
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
    
    if(verbose) message("--- ", prot_exp, " PROCESSED ")
    reproName <- paste0(output_name, ".qcplot.intensityStats.pdf")
    
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
  }# ends plotINTMISC
  
  if(plotPTMSTATS){
    if (prot_exp == "PH" | prot_exp == "UB") {
      if(verbose) message(">> GENERATING PTM ", prot_exp, " STATS ")
      modName <- paste0(output_name, "qcplot.ptmStats.pdf")
      
      x <-
        ggplot(evidencekeys, aes(x = BioReplicate, fill = MODIFICATION))
      x <-
        x + geom_bar(stat = "count",
                     position = position_dodge(width = 0.7),
                     width = 0.7)
      x <- x + theme_minimal()
      x <-
        x + theme(
          axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5
          ),
          legend.title = element_blank()
        )
      x <- x + ggtitle("Peptide Count in Biological Replicas")
      
      y <-
        ggplot(evidencekeys, aes(x = Condition, fill = MODIFICATION))
      y <-
        y + geom_bar(stat = "count",
                     position = position_dodge(width = 0.7),
                     width = 0.7)
      y <- y + theme_minimal()
      y <-
        y + theme(
          axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5
          ),
          legend.title = element_blank()
        )
      y <- y + ggtitle("Peptide Count in Conditions")
      
      u <-
        ggplot(evidencekeys,
               aes(x = BioReplicate, y = Intensity, fill = MODIFICATION))
      u <-
        u + geom_bar(stat = "identity",
                     position = position_dodge(width = 0.7),
                     width = 0.7)
      u <- u + theme_minimal()
      u <-
        u + theme(
          axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5
          ),
          legend.title = element_blank()
        )
      u <-
        u + ggtitle("Total Peptide Intensity in Biological Replicas")
      u <- u + scale_fill_brewer(palette = "Paired")
      
      z <-
        ggplot(evidencekeys,
               aes(x = Condition, y = Intensity, fill = MODIFICATION))
      z <-
        z + geom_bar(stat = "identity",
                     position = position_dodge(width = 0.7),
                     width = 0.7)
      z <- z + theme_minimal()
      z <-
        z + theme(
          axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5
          ),
          legend.title = element_blank()
        )
      z <- z + ggtitle("Total Peptide Intensity in Conditions")
      z <- z + scale_fill_brewer(palette = "Paired")
      
      #QC: Total Peptides
      if(printPDF) pdf(modName)
      print(x)
      print(y)
      print(u)
      print(z)
      if(printPDF) garbage <- dev.off()
    }
  } #plotPTMSTATS

  
  if(verbose) message(">> BASIC QUALITY CONTROL ANALYSIS COMPLETED!  ")
}
