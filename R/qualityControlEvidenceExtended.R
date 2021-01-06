
#' @title Extended Quality Control of the MaxQuant evidence.txt file
#'
#' @description Performs quality control based on the information available in
#' the MaxQuant `evidence.txt` file.
#' @param evidence_file (char or data.frame) The evidence file path and name, or
#' data.frame
#' @param keys_file (char or data.frame) The keys file path and name or
#' data.frame
#' @param output_name (char) prefix output name (no extension).
#' Default: "qcExtended_evidence"
#' @param isSILAC if `TRUE` processes SILAC input files. Default is `FALSE`
#' @param plotPSM (logical) `TRUE` generates peptide-spectrum-matches (PSMs) 
#' statistics plot: Page 1 shows the number of PSMs confidently identified 
#' in each BioReplicate. If replicates are present, Page 2 shows the mean 
#' number of PSMs per condition with error bar showing the standard error 
#' of the mean. Note that potential contaminant proteins are plotted separately. 
#' @param plotIONS (logical) `TRUE` generates peptide ion statistics plot: 
#' A peptide ion is defined in the context of m/z, in other words, an unique 
#' peptide sequence may give rise to multiple ions with different charge state 
#' and/or amino acid modification. Page 1 shows the number of ions confidently
#'  identified in each BioReplicate . If replicates are present, Page 2 shows 
#'  the mean number of peptide ions per condition with error bar showing the 
#'  standard error of the mean. Note that potential contaminant proteins are 
#'  plotted separately.  
#' @param plotTYPE (logical) `TRUE` generates identification type statistics 
#' plot: MaxQuant classifies each peptide identification into different 
#' categories (e.g., MSMS, MULTI-MSMS, MULTI-SECPEP). 
#' Page 1 shows the distribution of identification type in each BioReplicate
#' @param plotPEPTIDES (logical) `TRUE` generates peptide statistics plot: 
#' Page 1 shows the number of unique peptide sequences (disregard the charge 
#' state or amino acid modifications) confidently identified in each 
#' BioReplicate. If replicates are present, Page 2 shows the mean number of 
#' peptides per condition with error bar showing the standard error of the 
#' mean. Note that potential contaminant proteins are plotted separately. 
#' Pages 3 and 4 show peptide identification intersection between 
#' BioReplicates (the bars are ordered by degree or frequency, respectively),
#'  and Page 4 shows the intersections across conditions instead of 
#'  BioReplicates.
#' @param plotPEPTOVERLAP (logical) `TRUE` Show peptide identification 
#' intersection between BioReplicates and Conditions
#' @param plotPROTEINS (logical) `TRUE` generates protein statistics plot: 
#' Page 1 shows the number of protein groups confidently identified in each 
#' BioReplicate. If replicates are present, Page 2 shows the mean number of 
#' protein groups per condition with error bar showing the standard error of 
#' the mean. Note that potential contaminant proteins are plotted separately.
#' Pages 3 and 4 show peptide identification intersection between 
#' BioReplicates (the bars are ordered by degree or frequency, respectively), 
#' and Page 4 shows the intersections across conditions instead of 
#' BioReplicates.
#' @param plotPROTOVERLAP (logical) `TRUE` Show protein identification 
#' intersection between BioReplicates and Conditions
#' @param plotPIO (logical) `TRUE` generates oversampling statistics plot: 
#' Page 1 shows the proportion of all peptide ions (including peptides 
#' matched across runs) fragmented once, twice and thrice or more. 
#' Page 2 shows the  proportion of peptide ions (with intensity detected) 
#' fragmented once, twice and thrice or more. Page 3 shows the proportion of 
#' peptide ions (with intensity detected and MS/MS identification) fragmented 
#' once, twice and thrice or more
#' @param plotCS (logical) `TRUE` generates charge state plot: Page 1 shows 
#' the charge state distribution of PSMs confidently identified in each 
#' BioReplicate.
#' @param plotME (logical) `TRUE` generates precursor mass error plot: 
#' Page 1 shows the distribution of precursor error for all PSMs confidently 
#' identified in each BioReplicate.
#' @param plotMOCD (logical) `TRUE` generates precursor mass-over-charge plot: 
#' Page 1 shows the distribution of precursor mass-over-charge for all PSMs 
#' confidently identified in each BioReplicate. 
#' @param plotPEPICV (logical) `TRUE` generates peptide intensity coefficient 
#' of variance (CV) plot: The CV is calculated for each feature (peptide ion) 
#' identified in more than one replicate. Page 1 shows the distribution of 
#' CV's for each condition, while Page 2 shows the distribution of CV's 
#' within 4 bins of intensity (i.e., 4 quantiles of average intensity).
#' @param plotPEPDETECT (logical) `TRUE` generates peptide detection 
#' frequency plot: Page 1 summarizes the frequency that each peptide is 
#' detected across BioReplicates of each condition, showing the percentage 
#' of peptides detected once, twice, thrice, and so on (for whatever number 
#' of replicates each condition has).
#' @param plotPROTICV (logical) `TRUE` generates protein intensity coefficient 
#' of variance (CV) plot: The CV is calculated for each protein (after summing 
#' the peptide intensities) identified in more than one replicate. 
#' Page 1 shows the distribution of CV's for each condition, while Page 2 
#' shows the distribution of CV's within 4 bins of intensity (i.e., 4 
#' quantiles of average intensity).
#' @param plotPROTDETECT (logical) `TRUE` generates protein detection 
#' frequency plot: Page 1 summarizes the frequency that each protein group 
#' is detected across BioReplicates of each condition, showing the 
#' percentage of proteins detected once, twice, thrice, and so on (for 
#' whatever number of replicates each condition has). Page 2 shows the feature 
#' (peptide ion) intensity distribution within each BioReplicate (potential 
#' contaminant proteins are plot separately). Page 3 shows the density of 
#' feature intensity for different feature types (i.e., MULTI-MSMS, 
#' MULTI-SECPEP).
#' @param plotIDoverlap (logical) `TRUE` generates pairwise identification 
#' heatmap overlap: Pages 1 and 2 show pairwise peptide and protein overlap 
#' between any 2 BioReplicates, respectively. 
#' @param plotPCA (logical) `TRUE` generates PCA and pairwise intensity correlation: 
#' Page 1 and 3 show pairwise peptide and protein intensity correlation and 
#' scatter plot between any 2 BioReplicates, respectively. Page 2 and 4 show 
#' Principal Component Analysis at the intensity level for both peptide and 
#' proteins, respectively. 
#' @param plotSP (logical)  `TRUE` generates sample quality metrics: Page 1 
#' shows missing cleavage distribution of all peptides confidently identified 
#' in each BioReplicate. Page 2 shows the fraction of peptides with at least 
#' one methionine oxidized in each BioReplicate. 
#' @param printPDF If `TRUE` (default) prints out the pdfs. Warning: plot
#' objects are not returned due to the large number of them. 
#' @param verbose (logical) `TRUE` (default) shows function messages
#' @details all the plots are generated by default
#' @return A number of QC plots based on the evidence file
#' @keywords qc, evidence, keys
#' @examples
#' # Testing warning if files are not submitted
#' test <- artmsQualityControlEvidenceExtended(evidence_file = NULL,
#' keys_file = NULL)
#' @export
artmsQualityControlEvidenceExtended <- function(evidence_file,
                                                keys_file,
                                                output_name = "qcExtended_evidence",
                                                isSILAC = FALSE,
                                                plotPSM = TRUE,
                                                plotIONS = TRUE,
                                                plotTYPE = TRUE,
                                                plotPEPTIDES = TRUE,
                                                plotPEPTOVERLAP = TRUE,
                                                plotPROTEINS = TRUE,
                                                plotPROTOVERLAP = TRUE,
                                                plotPIO = TRUE,
                                                plotCS = TRUE,
                                                plotME = TRUE,
                                                plotMOCD = TRUE,
                                                plotPEPICV = TRUE,
                                                plotPEPDETECT = TRUE,
                                                plotPROTICV = TRUE,
                                                plotPROTDETECT = TRUE,
                                                plotIDoverlap = TRUE,
                                                plotPCA = TRUE,
                                                plotSP = TRUE,
                                                printPDF = TRUE,
                                                verbose = TRUE){
  
  # Local variables
  keysilac = NULL
  
  # DEBUG
  # evidence_file = artms_data_ph_evidence
  # keys_file = artms_data_ph_keys
  # output_name = "qcExtended_evidence"
  # isSILAC = FALSE
  # plotPSM = TRUE
  # plotIONS = TRUE
  # plotTYPE = TRUE
  # plotPEPTIDES = TRUE
  # plotPROTEINS = TRUE
  # plotPIO = TRUE
  # plotCS = TRUE
  # plotME = TRUE
  # plotMOCD = TRUE
  # plotPEPICV = TRUE
  # plotPEPDETECT = TRUE
  # plotPROTICV = TRUE
  # plotPROTDETECT = TRUE
  # plotIDoverlap = TRUE
  # plotPCA = TRUE
  # plotSP = TRUE
  # printPDF = FALSE
  # verbose = TRUE
  
  
  
  if(verbose){
    message("---------------------------------------------")
    message("artMS: EXTENDED QUALITY CONTROL (-evidence.txt based)")
    message("---------------------------------------------")
  }
  
  if (is.null(evidence_file) & is.null(keys_file)) {
    return("Evidence and keys cannot be NULL")
  }
  
  hmcol <- colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(256)
  
  # DATA PREPARATION----
  
  # OPEN KEYS
  keys <- .artms_checkIfFile(keys_file)
  keys <- .artms_checkRawFileColumnName(keys)
  
  # Check SILAC
  if(isSILAC){
    evidence_silac  <- artmsSILACtoLong(evidence_file,
                                        output = NULL,
                                        verbose = verbose)
    
    evidencekeys <- .artmsMergeSilacEvidenceKeys(evisilac = evidence_silac,
                                                 keysilac = keys)
  }else{
    evidencekeys <- artmsMergeEvidenceAndKeys(evidence_file, 
                                              keys,
                                              verbose = verbose)
  }
  
  colnames(evidencekeys) <- tolower(colnames(evidencekeys))
  
  if ( any(grepl("+", evidencekeys$reverse)) ) {
    evidencekeys <- subset(evidencekeys, reverse != "+")
  }
  
  # Add potential.contaminant column in case search was performed without it
  `%ni%` <- Negate(`%in%`)
  if ("potential.contaminant" %ni% colnames(evidencekeys)) {
    evidencekeys$potential.contaminant <- "no"
  }
  
  evidencekeys$potential.contaminant <- gsub("\\+", "yes", evidencekeys$potential.contaminant)
  evidencekeys$potential.contaminant <- gsub("^$", "no", evidencekeys$potential.contaminant)
  
  isFRACTION <- FALSE
  if("fraction" %in% colnames(evidencekeys)){
    if( length(unique(evidencekeys$fraction)) > 1 ){
      isFRACTION = TRUE
    }
  }
  
  evidencekeys.dt <- data.table::data.table(evidencekeys)
  
  evidence2 <- plyr::ddply(
    evidencekeys,
    c("bioreplicate", "condition", "potential.contaminant"),
    plyr::summarise,
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    PSMs = sum(ms.ms.count),
    Ions = length(sequence[ms.ms.count > 0]),
    Peptides = length(unique(sequence[ms.ms.count >
                                        0])),
    Proteins = length(unique(proteins[ms.ms.count >
                                        0])),
    Intensity = as.numeric(format(
      sum(intensity, na.rm = TRUE),
      scientific = TRUE,
      digits = 2
    )),
    PeakWidth.mean = mean(retention.length, na.rm =
                            TRUE)
  )
  
  if (isFRACTION) {
    evidence2fx <- plyr::ddply(
      evidencekeys,
      c(
        "bioreplicate",
        "condition",
        "potential.contaminant",
        "fraction"
      ),
      plyr::summarise,
      PSMs = sum(ms.ms.count),
      Ions = length(sequence[ms.ms.count > 0]),
      Peptides = length(unique(sequence[ms.ms.count >
                                          0])),
      Proteins = length(unique(proteins[ms.ms.count >
                                          0])),
      Intensity = as.numeric(format(
        sum(intensity, na.rm = TRUE),
        scientific = TRUE,
        digits = 2
      )),
      PeakWidth.mean = mean(retention.length, na.rm =
                              TRUE)
    )
  }
  
  evidence3 <- plyr::ddply(evidence2,
                           c("condition", "potential.contaminant"),
                           plyr::summarise,
                           PSMs.mean = mean(PSMs),
                           PSMs.max = max(PSMs),
                           PSMs.min = min(PSMs),
                           PSMs.sem = sd(PSMs) / sqrt(length(PSMs)),
                           Ions.mean = mean(Ions),
                           Ions.max = max(Ions),
                           Ions.min = min(Ions),
                           Ions.sem = sd(Ions) / sqrt(length(Ions)),
                           Peptides.mean = mean(Peptides),
                           Peptides.max = max(Peptides),
                           Peptides.min = min(Peptides),
                           Peptides.sem = sd(Peptides) / sqrt(length(Peptides)),
                           Proteins.mean = mean(Proteins),
                           Proteins.max = max(Proteins),
                           Proteins.min = min(Proteins),
                           Proteins.sem = sd(Proteins) / sqrt(length(Proteins)),
                           Intensity.mean = mean(Intensity),
                           Intensity.max = max(Intensity),
                           Intensity.min = min(Intensity),
                           Intensity.sem = sd(Intensity) / sqrt(length(Intensity)),
                           PeakWidth.mean.mean = mean(PeakWidth.mean),
                           PeakWidth.mean.max = max(PeakWidth.mean),
                           PeakWidth.mean.min = min(PeakWidth.mean),
                           PeakWidth.mean.sem = sd(PeakWidth.mean) / sqrt(length(PeakWidth.mean))
                           )
    
  
  # - Oversampling
  ## Oversampling is calculated for peptides identified by MS/MS and detected
  ## in MS1 (i.e., peptides with intensity calculated)
  oversampling0 <- evidencekeys.dt[, .N, by = list(ms.ms.count, bioreplicate)]
  oversampling0$MSMS.counts <- findInterval(oversampling0$ms.ms.count, seq(1, 3, by = 1))
  oversampling0 <- oversampling0[, list(N = sum(N)), by = list(MSMS.counts, bioreplicate)]
  oversampling0$MSMS.counts <- paste("n=", oversampling0$MSMS.counts, sep = "")
  oversampling0$MSMS.counts <- gsub("n=3", "n=3+", oversampling0$MSMS.counts)
  oversampling0.total <- evidencekeys.dt[, .N, by = list(bioreplicate)]
  oversampling0 <- merge(oversampling0,
                         oversampling0.total,
                         by = "bioreplicate",
                         all = TRUE)
    
  oversampling0$FxOverSamp <- as.numeric(format(100 * (oversampling0$N.x / oversampling0$N.y), digits = 3))
  
  oversampling1 <- evidencekeys.dt[intensity > 0, .N, by = list(ms.ms.count, bioreplicate)]
  oversampling1$MSMS.counts <- findInterval(oversampling1$ms.ms.count, seq(1, 3, by = 1))
  oversampling1 <- oversampling1[, list(N = sum(N)), by = list(MSMS.counts, bioreplicate)]
  oversampling1$MSMS.counts <- paste("n=", oversampling1$MSMS.counts, sep = "")
  oversampling1$MSMS.counts <- gsub("n=3", "n=3+", oversampling1$MSMS.counts)
  oversampling1.total <- evidencekeys.dt[intensity > 0, .N, by = list(bioreplicate)]
  oversampling1 <- merge(oversampling1,
                         oversampling1.total,
                         by = "bioreplicate",
                         all = TRUE)
    
  oversampling1$FxOverSamp <- as.numeric(format(100 * (oversampling1$N.x / oversampling1$N.y), digits = 3))
  
  oversampling2 <- evidencekeys.dt[intensity > 0 &
                                     ms.ms.count > 0, .N, by = list(ms.ms.count, bioreplicate)]
    
  oversampling2$MSMS.counts <- findInterval(oversampling2$ms.ms.count, seq(1, 3, by = 1))
  oversampling2 <- oversampling2[, list(N = sum(N)), by = list(MSMS.counts, bioreplicate)]
  oversampling2$MSMS.counts <- paste("n=", oversampling2$MSMS.counts, sep = "")
  oversampling2$MSMS.counts <- gsub("n=3", "n=3+", oversampling2$MSMS.counts)
  oversampling2.total <- evidencekeys.dt[intensity > 0 & ms.ms.count > 0, .N, by = list(bioreplicate)]
  oversampling2 <- merge(oversampling2,
                         oversampling2.total,
                         by = "bioreplicate",
                         all = TRUE)
  oversampling2$FxOverSamp <- as.numeric(format(100 * (oversampling2$N.x / oversampling2$N.y), digits = 3))
  
  # Charge state
  chargeState <- evidencekeys.dt[, .N, by = list(charge, bioreplicate)]
  chargeState$charge <- paste("z=", chargeState$charge, sep = "")
  chargeState.total <- evidencekeys.dt[, .N, by = list(bioreplicate)]
  chargeState <- merge(chargeState,
                       chargeState.total,
                       by = "bioreplicate",
                       all = TRUE)
    
  chargeState$FxOverSamp <- as.numeric(format(100 * (chargeState$N.x / chargeState$N.y), digits = 1))
  
  chargeStateCond <- evidencekeys.dt[, .N, by = list(charge, condition)]
  chargeStateCond$charge <- paste("z=", chargeStateCond$charge, sep = "")
  chargeStateCond.total <- evidencekeys.dt[, .N, by = list(condition)]
  chargeStateCond <- merge(chargeStateCond,
                           chargeStateCond.total,
                           by = "condition",
                           all = TRUE)
  chargeStateCond$FxOverSamp <- as.numeric(format(100 * (chargeStateCond$N.x / chargeStateCond$N.y), digits = 1))
    
  
  # Type
  mstype <- evidencekeys.dt[, .N, by = list(type, bioreplicate, condition)]
  mstype.total <- evidencekeys.dt[, .N, by = list(bioreplicate)]
  mstype <- merge(mstype, mstype.total, by = "bioreplicate", all = TRUE)
  mstype$FxOverSamp <- as.numeric(format(100 * (mstype$N.x / mstype$N.y), digits = 2))
  
  # PLOTS-----
  if(verbose) message(">> GENERATING QC PLOTS ")
  
  nsamples <- length(unique(evidencekeys$bioreplicate))
  nconditions <- length(unique(evidencekeys$condition))
  
  ### PSM
  if (plotPSM) {
    if(verbose) message("--- Plot PSM", appendLF = FALSE)
    if(printPDF) pdf(
      paste0(output_name,'.qcplot.PSM.pdf'),
      width = 10, #nsamples * 3
      height = 6,
      onefile = TRUE
    )
    
    aa <- ggplot(evidence2,
                 aes(
                   x = bioreplicate,
                   y = PSMs,
                   fill = factor(condition)
                 )) +
      geom_bar(stat = "identity",
               position = position_dodge(width = 1),
               alpha = 0.7, na.rm = TRUE) +
      geom_text(aes(label = round(PSMs, digits = 0)),
                hjust = 0.5,
                vjust = -0.5,
                size = 2,
                angle = 0,
                position = position_dodge(width = 1)
                ) +
      facet_wrap( ~ potential.contaminant, ncol = 1) +
      xlab("Experiment") + ylab("Counts") +
      labs(title = "Number of Peptide-spectrum-matches (PSMs)",           
           caption = "bottom = Potential contaminants; top = non-contaminants",
           fill = "") +
      theme_linedraw() +
      theme(legend.text = element_text(size = 8)) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 0,
        size = 8
      )) +
      theme(axis.text.y = element_text(size = 10)) +
      theme(axis.title.y = element_text(size = 10)) +
      theme(plot.title = element_text(size = 12)) +
      
      scale_fill_brewer(palette = "Spectral")
    
    print(aa)
    
    ab <- ggplot(evidence3,
                 aes(
                   x = condition,
                   y = PSMs.mean,
                   fill = factor(potential.contaminant)
                 )) +
      geom_bar(stat = "identity",
               position = position_dodge(width = 1),
               alpha = 0.7, na.rm = TRUE) +
      geom_errorbar(
        aes(ymin = PSMs.mean - PSMs.sem, ymax = PSMs.mean + PSMs.sem),
        width = .2,
        position = position_dodge(.9)
      ) +
      geom_text(
        aes(label = round(PSMs.mean, digits = 0)),
        hjust = 0.5,
        vjust = -0.5,
        size = 2,
        angle = 0,
        position = position_dodge(width = 1)
      ) +
      xlab("Condition") + ylab("Counts") +
      theme_linedraw() +
      labs(title = "Mean number of PSMs per condition",
           caption = "yes = Potential contaminants; no = non-contaminants",
           fill = "Conditions") +
      theme(legend.text = element_text(size = 8)) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 8
      )) +
      theme(axis.text.y = element_text(size = 10)) +
      theme(axis.title.y = element_text(size = 10)) +
      theme(plot.title = element_text(size = 12)) +
      scale_fill_brewer(palette = "Spectral")
    
    print(ab)
    
    if (isFRACTION) {
      ac <- ggplot(evidence2fx, aes(
        x = bioreplicate,
        y = PSMs,
        fill = factor(fraction)
      )) +
        geom_bar(stat = "identity", alpha = 0.7, na.rm = TRUE) +
        geom_text(
          aes(label = round(PSMs, digits = 0)),
          # hjust = 0.5,
          # vjust = 1.5,
          # size = 2,
          # position = position_stack()
          hjust = 1,
          size = 2.7,
          angle = 90,
          position = position_dodge(width = 1)
        ) +
        facet_wrap( ~ potential.contaminant,
                    ncol = 1,
                    scales = "free") +
        xlab("Experiment") + ylab("Counts") +
        labs(title = "Number of PSMs per Fraction",
             caption = "bottom = Potential contaminants; top = non-contaminants",
             fill = "") +
        theme(legend.text = element_text(size = 10)) +
        theme(axis.text.x = element_text(
          angle = 90,
          hjust = 0,
          size = 10
        )) +
        theme(axis.text.y = element_text(size = 10)) +
        theme(axis.title.y = element_text(size = 10)) +
        theme(plot.title = element_text(size = 12))
        
      print(ac)
    }
    
    if(printPDF) garbage <- dev.off()
    if(verbose) message(" done ")
  }
  
  
  ### IONS
  if (plotIONS) {
    if(verbose) message("--- Plot IONS", appendLF = FALSE)
    
    if(printPDF) pdf(
      paste0(output_name,'.qcplot.Ions.pdf'),
      width = 10, #nsamples * 3
      height = 6,
      onefile = TRUE
    )
    
    ba <- ggplot(evidence2,
                 aes(
                   x = bioreplicate,
                   y = Ions,
                   fill = factor(condition)
                 )) +
      geom_bar(stat = "identity",
               position = position_dodge(width = 1),
               alpha = 0.7, na.rm = TRUE) +
      geom_text(
        aes(label = round(Ions, digits = 0)),
        hjust = 0.5,
        vjust = -0.5,
        size = 2,
        position = position_dodge(width = 1)
      ) +
      facet_wrap( ~ potential.contaminant, ncol = 1) +
      xlab("Experiment") + ylab("Counts") +
      labs(title = "Number of unique Peptide Ions",
           caption = "bottom = Potential contaminants; top = non-contaminants",
           fill = "") +
      theme_linedraw() +
      theme(legend.text = element_text(size = 8)) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 8
      )) +
      theme(axis.text.y = element_text(size = 10)) +
      theme(axis.title.y = element_text(size = 10)) +
      theme(plot.title = element_text(size = 12)) +
      scale_fill_brewer(palette = "Spectral")
    print(ba)
    
    bb <- ggplot(evidence3, aes(
      x = condition,
      y = Ions.mean,
      fill = factor(potential.contaminant)
    )) +
      geom_bar(stat = "identity",
               position = position_dodge(width = 1),
               alpha = 0.7, na.rm = TRUE) +
      geom_errorbar(
        aes(ymin = Ions.mean - Ions.sem, ymax = Ions.mean + Ions.sem),
        width = .2,
        position = position_dodge(.9)
      ) +
      geom_text(
        aes(label = round(Ions.mean, digits = 0)),
        hjust = 0.5,
        vjust = -0.5,
        size = 2,
        position = position_dodge(width = 1)
      ) +
      xlab("Condition") + ylab("Counts") +
      labs(title = "Mean number of unique Peptide Ions",
           subtitle = "Error bar = std error of the mean",
           fill = "Contaminants") +
      theme_linedraw() +
      theme(legend.text = element_text(size = 8)) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 8
      )) +
      theme(axis.text.y = element_text(size = 10)) +
      theme(axis.title.y = element_text(size = 10)) +
      theme(plot.title = element_text(size = 12)) +
      scale_fill_brewer(palette = "Set1")
      
    print(bb)
    
    if (isFRACTION) {
      bc <- ggplot(evidence2fx, aes(
        x = bioreplicate,
        y = Ions,
        fill = factor(fraction)
      )) +
        geom_bar(stat = "identity", alpha = 0.7, na.rm = TRUE) +
        geom_text(
          aes(label = round(Ions, digits = 0)),
          hjust = 0.5,
          vjust = 1.5,
          size = 2,
          angle = 90,
          position = position_dodge(width = 1)
        ) +
        facet_wrap( ~ potential.contaminant,
                    ncol = 1,
                    scales = "free") +
        xlab("Experiment") + ylab("Counts") +
        labs(title = "Number of unique Peptide Ions in each Fraction",
             caption = "bottom = Potential contaminants; top = non-contaminants",
             fill = "Fractions") +
        theme_linedraw() +
        theme(legend.text = element_text(size = 10)) +
        theme(axis.text.x = element_text(
          angle = 90,
          hjust = 0,
          size = 8
        )) +
        theme(axis.text.y = element_text(size = 10)) +
        theme(axis.title.y = element_text(size = 10)) +
        theme(plot.title = element_text(size = 12))
        
      print(bc)
    }
    
    if(printPDF) garbage <- dev.off()
    
    if(verbose) message(" done ")
  }
  
  #### TYPE ######
  if (plotTYPE) {
    
    if(verbose) message("--- Plot TYPE", appendLF = FALSE)
    
    if(printPDF) pdf(
      paste0(output_name,'.qcplot.Type.pdf'),
      width = 10, #nsamples * 3
      height = 6,
      onefile = TRUE
    )
    
    ca <- ggplot(mstype, aes(
      x = bioreplicate,
      y = FxOverSamp,
      fill = factor(type))) +
      geom_bar(stat = "identity",
               position = position_stack(),
               alpha = 0.7, na.rm = TRUE) +
      ggrepel::geom_text_repel(
        aes(label = round(FxOverSamp, digits = 1)),
        vjust = 1,
        size = 2,
        position = position_stack()
      ) +
      xlab("Experiment") + ylab("Fraction") +
      labs(title = "Type of MaxQuant Feature Identification",
           fill = "Type") +
      theme_linedraw() +
      theme(legend.text = element_text(size = 8)) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 8
      )) +
      theme(axis.text.y = element_text(size = 10)) +
      theme(axis.title.y = element_text(size = 10)) +
      theme(plot.title = element_text(size = 12)) +
      scale_fill_brewer(palette = "Set1")
      
    print(ca)
    
    if(printPDF) garbage <- dev.off()
    
    if(verbose) message(" done ")
  }
  
  ### PEPTIDES ###
  if (plotPEPTIDES) {
    if(verbose) message("--- Plot PEPTIDES", appendLF = FALSE)
    
    if(printPDF) pdf(
      paste0(output_name,'.qcplot.Peptides.pdf'),
      width = 10, #nsamples * 3
      height = 6,
      onefile = TRUE
    )
    
    da <- ggplot(evidence2, aes(
      x = bioreplicate,
      y = Peptides,
      fill = factor(condition))) +
      geom_bar(stat = "identity",
               position = position_dodge(width = 1),
               alpha = 0.7, na.rm = TRUE) +
      geom_text(
        aes(label = round(Peptides, digits = 0)),
        hjust = 0.5,
        vjust = -0.5,
        size = 2,
        angle = 0,
        position = position_dodge(width = 1)
      ) +
      facet_wrap( ~ potential.contaminant, ncol = 1) +
      xlab("Experiment") + ylab("Counts") +
      labs(title = "Number of unique Peptides",
           caption = "bottom = Potential contaminants; top = non-contaminants",
           fill = "Condition") +
      theme_linedraw() +
      theme(legend.text = element_text(size = 8)) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 8
      )) +
      theme(axis.text.y = element_text(size = 10)) +
      theme(axis.title.y = element_text(size = 10)) +
      theme(plot.title = element_text(size = 12)) +
      scale_fill_brewer(palette = "Spectral")
      
    print(da)
    
    db <- ggplot(evidence3,
                 aes(
                   x = condition,
                   y = Peptides.mean,
                   fill = factor(potential.contaminant)
                 )) +
      geom_bar(stat = "identity",
               position = position_dodge(width = 1),
               alpha = 0.7, na.rm = TRUE) +
      geom_errorbar(
        aes(
          ymin = Peptides.mean - Peptides.sem,
          ymax = Peptides.mean + Peptides.sem
        ),
        width = .2,
        position = position_dodge(.9)
      ) +
      geom_text(
        aes(label = round(Peptides.mean, digits = 0)),
        hjust = 0.5,
        vjust = -0.5,
        size = 2,
        position = position_dodge(width = 1)
      ) +
      xlab("Condition") + ylab("Counts") +
      labs(title = "Number of unique Peptides",
           fill = "Contaminants") +
      theme_linedraw() +
      theme(legend.text = element_text(size = 8)) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 8
      )) +
      theme(axis.text.y = element_text(size = 10)) +
      theme(axis.title.y = element_text(size = 10)) +
      theme(plot.title = element_text(size = 12)) +
      scale_fill_brewer(palette = "Set1")
      
    print(db)
    
    if (isFRACTION) {
      dc <- ggplot(evidence2fx,
                   aes(
                     x = bioreplicate,
                     y = Peptides,
                     fill = factor(fraction)
                   )) +
        geom_bar(stat = "identity", alpha = 0.7, na.rm = TRUE) +
        geom_text(
          aes(label = round(Peptides, digits = 0)),
          # hjust = 0.5,
          # vjust = 1.5,
          # size = 7,
          # position = position_stack()
          hjust = 1,
          size = 2.7,
          angle = 90,
          position = position_dodge(width = 1)
        ) +
        facet_wrap( ~ potential.contaminant,
                    ncol = 1,
                    scales = "free") +
        xlab("Experiment") + ylab("Counts") +
        labs(title = "Number of unique Peptides",
             fill = "Contaminants") +
        theme(legend.text = element_text(size = 10)) +
        theme(axis.text.x = element_text(
          angle = 90,
          hjust = 0,
          size = 10
        )) +
        theme(axis.text.y = element_text(size = 10)) +
        theme(axis.title.y = element_text(size = 10)) +
        theme(plot.title = element_text(size = 12)) +
        theme_linedraw()
        
      print(dc)
    }
    
    if(printPDF) garbage <- dev.off()
    if(verbose) message(" done ")
  }
  
  ### PEPTIDES OVERLAP
  if (plotPEPTOVERLAP) {
    if(verbose) message("--- Plot PEPTIDE OVERLAP", appendLF = FALSE)
    
    lst.pepExp <- list()
    for (i in unique(evidencekeys$bioreplicate)) {
      lst.pepExp[[i]] <-
        unique(subset(evidencekeys, ms.ms.count > 0 &
                        bioreplicate == i)[, c("sequence")])
    }
    
    lst.pepCond <- list()
    for (i in unique(evidencekeys$condition)) {
      lst.pepCond[[i]] <- unique(subset(evidencekeys,
                                        ms.ms.count > 0 &
                                          condition == i)[, c("sequence")])
    }
    
    uspep1 <- UpSetR::upset(
      fromList(lst.pepExp),
      nsets = length(unique(evidencekeys$bioreplicate)),
      nintersects = 50,
      order.by = "freq",
      text.scale = c(1, 1, 1, 1, 0.75, 0.75)
    )
    
    uspep2 <- UpSetR::upset(
      fromList(lst.pepCond),
      nsets = length(unique(evidencekeys$condition)),
      nintersects = 50,
      order.by = "freq",
      text.scale = c(1, 1, 1, 1, 0.75, 0.75)
    )
    
    if(printPDF) pdf(
      paste0(output_name,'.qcplot.PeptidesOverlap.pdf'),
      width = 10,
      height = 6,
      onefile = TRUE
    )
    
    print(uspep1)
    print(uspep2)
    
    if(printPDF) garbage <- dev.off()
    if(verbose) message(" done ")
  }
  
  
  ## PROTEINS
  if (plotPROTEINS) {
    if(verbose) message("--- Plot PROTEINS", appendLF = FALSE)
    
    if(printPDF) pdf(
      paste0(output_name,'.qcplot.Proteins.pdf'),
      width = 10, #nsamples * 3
      height = 6,
      onefile = TRUE
    )
    
    ea <- ggplot(evidence2, aes(
      x = bioreplicate,
      y = Proteins,
      fill = factor(condition))) +
      geom_bar(stat = "identity",
               position = position_dodge(width = 1),
               alpha = 0.7, na.rm = TRUE) +
      geom_text(
        aes(label = round(Proteins, digits = 0)),
        hjust = 0.5,
        vjust = -0.5,
        size = 2,
        position = position_dodge(width = 1),
        angle = 0
      ) +
      facet_wrap( ~ potential.contaminant, ncol = 1) +
      xlab("Experiment") + ylab("Counts") +
      labs(title = "Number of unique Protein Groups",
           fill = "") +
      theme_linedraw() +
      theme(legend.text = element_text(size = 8)) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 8
      )) +
      theme(axis.text.y = element_text(size = 10)) +
      theme(axis.title.y = element_text(size = 10)) +
      theme(plot.title = element_text(size = 12)) +
      scale_fill_brewer(palette = "Spectral")
      
    print(ea)
    
    eb <- ggplot(evidence3,
                 aes(x = condition,
                     y = Proteins.mean,
                     fill = factor(potential.contaminant))) +
      geom_bar(stat = "identity",
               position = position_dodge(width = 1),
               alpha = 0.7, na.rm = TRUE) +
      geom_errorbar(
        aes(
          ymin = Proteins.mean - Proteins.sem,
          ymax = Proteins.mean + Proteins.sem
        ),
        width = .2,
        position = position_dodge(.9)
      ) +
      geom_text(
        aes(label = round(Proteins.mean, digits = 0)),
        hjust = 0.5,
        vjust = -0.5,
        size = 2,
        position = position_dodge(width = 1)
      ) +
      xlab("Condition") + ylab("Counts") +
      labs(title = "Mean number of unique Proteins",
           subtitle = "error bar = std error of the mean",
           fill = "Contaminants") +
      theme_linedraw() +
      theme(legend.text = element_text(size = 8)) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 8
      )) +
      theme(axis.text.y = element_text(size = 10)) +
      theme(axis.title.y = element_text(size = 10)) +
      theme(plot.title = element_text(size = 12)) +
      scale_fill_brewer(palette = "Set1")
      
    print(eb)
    
    if(printPDF) garbage <- dev.off()
    
    if(verbose) message(" done ")
  }

  ### PEPTIDES OVERLAP
  if (plotPROTOVERLAP) {
    if(verbose) message("--- Plot PROTEIN OVERLAP", appendLF = FALSE)
    
    lst.protExp <- list()
    for (i in unique(evidencekeys$bioreplicate)) {
      lst.protExp[[i]] <-
        unique(subset(evidencekeys, ms.ms.count > 0 &
                        bioreplicate == i)[, c("proteins")])
    }
    
    lst.protCond <- list()
    for (i in unique(evidencekeys$condition)) {
      lst.protCond[[i]] <-
        unique(subset(evidencekeys, ms.ms.count > 0 &
                        condition == i)[, c("proteins")])
    }
    
    upprot1 <- UpSetR::upset(
      fromList(lst.protExp),
      nsets = length(unique(evidencekeys$bioreplicate)),
      nintersects = 50,
      order.by = "freq",
      text.scale = c(1, 1, 1, 1, 0.75, 0.75)
    )
    
    upprot2 <- UpSetR::upset(
      fromList(lst.protCond),
      nsets = length(unique(evidencekeys$condition)),
      nintersects = 50,
      order.by = "freq",
      text.scale = c(1, 1, 1, 1, 0.75, 0.75)
    )
    
    if(printPDF) pdf(
      paste0(output_name,'.qcplot.ProteinOverlap.pdf'),
      width = 10,
      height = 6,
      onefile = TRUE
    )
    
    print(upprot1)
    print(upprot2)
    
    if(printPDF) garbage <- dev.off()
    if(verbose) message(" done ")
  }
  
  
  # Peptide ion oversampling
  if (plotPIO) {
    if(verbose) message("--- Plot Plot Ion Oversampling", appendLF = FALSE)
    
    if(printPDF) pdf(
      paste0(output_name,'.qcplot.PepIonOversampling.pdf'),
      width = 10, #nsamples * 3
      height = 6,
      onefile = TRUE
    )
    fa <- ggplot(
      oversampling0[with(oversampling0, order(MSMS.counts)), ],
      aes(
        x = bioreplicate,
        y = FxOverSamp,
        fill = MSMS.counts,
        label = FxOverSamp
      )
    ) +
      geom_bar(stat = "identity", alpha = 0.7, na.rm = TRUE) +
      geom_text(
        aes(label = round(FxOverSamp, digits = 1)),
        vjust = 1 ,
        size = 2,
        position = position_stack()
      ) +
      xlab("Experiment") + ylab("Fraction (percentage)") +
      labs(title = "Peptide ion oversampling",
           subtitle = "Based on all the peptides reported by MaxQuant") +
      theme_linedraw() +
      theme(legend.text = element_text(size = 8)) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 8
      )) +
      theme(axis.text.y = element_text(size = 10)) +
      theme(axis.title.y = element_text(size = 10)) +
      theme(plot.title = element_text(size = 12)) +
      scale_fill_brewer(palette = "Set1")
      
    print(fa)
    
    fb <- ggplot(oversampling1[with(oversampling1, order(MSMS.counts)), ],
                 aes(x = bioreplicate,
                     y = FxOverSamp,
                     fill = MSMS.counts,
                     label = FxOverSamp)) +
      geom_bar(stat = "identity", alpha = 0.7, na.rm = TRUE) +
      geom_text(
        aes(label = round(FxOverSamp, digits = 1)),
        vjust = 1 ,
        size = 2,
        position = position_stack()
      ) +
      xlab("Experiment") + ylab("Fraction (percentage)") +
      labs(title = "Peptide ion oversampling",
           subtitle = "Only peptides detected (MS1 AUC calculated)") +
      theme_linedraw() +
      theme(legend.text = element_text(size = 8)) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 8
      )) +
      theme(axis.text.y = element_text(size = 10)) +
      theme(axis.title.y = element_text(size = 10)) +
      theme(plot.title = element_text(size = 12)) +
      scale_fill_brewer(palette = "Set1")
      
    print(fb)
    
    fc <- ggplot(oversampling2[with(oversampling2, order(MSMS.counts)), ],
                 aes(
                   x = bioreplicate,
                   y = FxOverSamp,
                   fill = MSMS.counts,
                   label = FxOverSamp)) +
      geom_bar(stat = "identity", alpha = 0.7, na.rm = TRUE) +
      geom_text(
        aes(label = round(FxOverSamp, digits = 1)),
        vjust = 1 ,
        size = 2,
        position = position_stack()
      ) +
      xlab("Experiment") + ylab("Fraction (percentage)") +
      labs(title = "Peptide ion oversampling",
           subtitle = "Only peptides detected (MS1 AUC calculated) and identified (confidence MS/MS)") +
      theme_linedraw() +
      theme(legend.text = element_text(size = 8)) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 8
      )) +
      theme(axis.text.y = element_text(size = 10)) +
      theme(axis.title.y = element_text(size = 10)) +
      theme(plot.title = element_text(size = 12)) +
      scale_fill_brewer(palette = "Set1")
      
    print(fc)
    
    if(printPDF) garbage <- dev.off()
    
    if(verbose) message(" done ")
  }
  
  
  # Charge state distribution
  if (plotCS) {
    if(verbose) message("--- Plot Charge State", appendLF = FALSE)
    
    if(printPDF) pdf(
      paste0(output_name,'.qcplot.ChargeState.pdf'),
      width = 10, #nsamples * 3
      height = 6,
      onefile = TRUE)
    
    ga <- ggplot(chargeState[with(chargeState, order(charge)), ],
                 aes(x = bioreplicate, y = FxOverSamp, fill = charge)) +
      geom_bar(stat = "identity", alpha = 0.7, na.rm = TRUE) +
      theme_linedraw() +
      ggrepel::geom_text_repel(
        aes(label = round(FxOverSamp, digits = 1)),
        vjust = -0.5,
        size = 2,
        position = position_stack()
      ) +
      # geom_text(
      #   aes(label = round(FxOverSamp, digits = 1)),
      #   vjust = 1 ,
      #   size = 2,
      #   position = position_stack()
      # ) +
      xlab("Experiment") + ylab("Fraction (percentage)") +
      ggtitle("Precursor charge state distribution") +
      theme_linedraw() +
      theme(legend.text = element_text(size = 8)) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 8
      )) +
      theme(axis.text.y = element_text(size = 10)) +
      theme(axis.title.y = element_text(size = 10)) +
      theme(plot.title = element_text(size = 12)) +
      scale_fill_brewer(palette = "Spectral")
    
    print(ga)
    
    if(printPDF) garbage <- dev.off()
    
    if(verbose) message(" done ")
  }
  
  # Mass error
  if (plotME) {
    if(verbose) message("--- Plot Mass Error", appendLF = FALSE)
    if(printPDF) pdf(
      paste0(output_name,'.qcplot.MassError.pdf'),
      width = 10, #nsamples * 3
      height = 6,
      onefile = TRUE
    )
    fa <- ggplot(evidencekeys,
                 aes(x = bioreplicate, y = uncalibrated.mass.error..ppm.)) +
      geom_boxplot(varwidth = TRUE, 
                   aes(fill = factor(condition)), 
                   alpha = 0.7,
                   na.rm = TRUE) +
      geom_text(
        data = aggregate(
          uncalibrated.mass.error..ppm. ~ bioreplicate,
          evidencekeys,
          median
        ),
        aes(
          label = round(uncalibrated.mass.error..ppm., digits = 1),
          y = max(evidencekeys$uncalibrated.mass.error..ppm., na.rm = TRUE) + 2
        ),
        size = 2
      ) +
      xlab("Experiment") + ylab("mass error") +
      labs(title = "Precursor mass error (in ppm) distribution",
           caption = "Global median mass error on the top",
           fill = "") +
      theme_linedraw() +
      theme(legend.text = element_text(size = 8)) +
      theme(axis.text.x = element_text(angle = 90,
                                       hjust = 1,
                                       size = 8)) +
      theme(axis.text.y = element_text(size = 10)) +
      theme(axis.title.y = element_text(size = 10)) +
      theme(plot.title = element_text(size = 12)) +
      scale_fill_brewer(palette = "Spectral")
    print(fa)
    if(printPDF) garbage <- dev.off()
    if(verbose) message(" done ")
  }
  
  # mass-over-charge distribution
  if (plotMOCD) {
    if(verbose) message("--- Plot Mass-over-Charge distribution", appendLF = FALSE)
    
    if(printPDF) pdf(
      paste0(output_name,'.qcplot.MZ.pdf'),
      width = 10, #nsamples * 3
      height = 6,
      onefile = TRUE
    )
    ga <- ggplot(evidencekeys, aes(x = bioreplicate, y = m.z)) +
      geom_boxplot(varwidth = TRUE, aes(fill = factor(condition)), alpha = 0.7, na.rm = TRUE) +
      geom_text(
        data = aggregate(m.z ~ bioreplicate, evidencekeys, median),
        aes(label = round(m.z, digits = 1),
            y = max(evidencekeys$m.z, na.rm = TRUE) + 30),
        size = 2,
        angle = 0) +
      xlab("Experiment") + ylab("m/z") +
      labs(title = "Precursor mass-over-charge distribution",
           caption = "Global median m/z on the top",
           fill = "") +
      theme_linedraw() +
      theme(legend.text = element_text(size = 8)) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 8
      )) +
      theme(axis.text.y = element_text(size = 10)) +
      theme(axis.title.y = element_text(size = 10)) +
      theme(plot.title = element_text(size = 12)) +
      scale_fill_brewer(palette = "Spectral")
    print(ga)
    
    if(printPDF) garbage <- dev.off()
    if(verbose) message(" done ")
  }
  
  
  #m Peptide Intensity CV
  if (plotPEPICV) {
    if(verbose) message("--- Plot Peptide Intensity CV", appendLF = FALSE)
    if(printPDF) pdf(
      paste0(output_name,'.qcplot.PeptideIntensityCV.pdf'),
      width = 10, #nsamples * 3
      height = 6,
      onefile = TRUE
    )
    peptCV <- data.table::data.table(subset(evidencekeys,!is.na(intensity)))
    peptCV$intensity <- as.numeric(peptCV$intensity)
  
    peptCV <- peptCV[, list(intensity = sum(intensity, na.rm = TRUE)),
                     by = list(condition, bioreplicate, modified.sequence)]
      
    peptCV <- peptCV[, list(pCV = 100 * (sd(intensity) / mean(intensity)),
                            sumInt = sum(intensity),
                            pDX = length(unique(bioreplicate))),
                     by = list(condition, modified.sequence)]
    peptCV <- peptCV[, bin.all := .artms_qcut(sumInt, 4)]
    peptCV <- peptCV[, bin.condition := .artms_qcut(sumInt, 4),
                     by = condition]
      
    
    ha <- ggplot(subset(peptCV,!is.na(pCV)), aes(condition, pCV)) +
      geom_boxplot(varwidth = TRUE, aes(fill = factor(condition)), alpha = 0.7, na.rm = TRUE) +
      geom_text(
        data = aggregate(pCV ~ condition, subset(peptCV,!is.na(pCV)), median),
        aes(label = round(pCV, digits = 1),
            y = max(peptCV$pCV, na.rm = TRUE) + 1),
        size = 2.5) +
      geom_text(
        data = aggregate(pCV ~ condition, subset(peptCV,!is.na(pCV)), length),
        aes(label = round(pCV, digits = 1), y = 0),
        size = 2) +
      xlab("Condition") + ylab("Coefficient of variance (%)") +
      labs(title = "Distribution of peptide feature intensity CV",
           fill = "") +
      theme_linedraw() +
      theme(legend.text = element_text(size = 8)) +
      theme(axis.text.x = element_text(angle = 90, 
                                       size = 8)) +
      theme(axis.text.y = element_text(size = 10)) +
      theme(axis.title.y = element_text(size = 10)) +
      theme(plot.title = element_text(size = 12)) +
      scale_fill_brewer(palette = "Spectral")
      
    print(ha)
    
    hb <- ggplot(subset(peptCV,!is.na(pCV)), aes(interaction(condition, bin.condition), pCV)) +
      geom_boxplot(varwidth = TRUE, aes(fill = factor(condition)), alpha = 0.7, na.rm = TRUE) +
      geom_text(
        data = aggregate(
          pCV ~ condition + bin.condition,
          subset(peptCV,!is.na(pCV)),
          median
        ),
        aes(
          label = round(pCV, digits = 1),
          y = max(peptCV$pCV, na.rm = TRUE) + 1
        ),
        size = 2,
        angle = 0
      ) +
      geom_text(
        data = aggregate(
          pCV ~ condition + bin.condition,
          subset(peptCV,!is.na(pCV)),
          length
        ),
        aes(label = round(pCV, digits = 1), y = 0),
        size = 2.5,
        angle = 0
      ) +
      xlab("Condition") + ylab("Coefficient of variance (%)") +
      labs(title = "Distribution of peptide feature intensity CV",
           caption = "For each condition, peptides were ranked by summed intensity 
           and the CV for each peptide was calculated. 
           Each condition shows 4 distribution (box) for low (1) to high (4) intensity peptides. 
           Overall median CV within each bin/condition is shown on the top and 
           number of features used to calculate CV is given on the bottom",
           fill = "") +
      theme_linedraw() +
      theme(legend.text = element_text(size = 8)) +
      theme(axis.text.x = element_text(angle = 90, 
                                       size = 8)) +
      theme(axis.text.y = element_text(size = 10)) +
      theme(axis.title.y = element_text(size = 10)) +
      theme(plot.title = element_text(size = 12, hjust = 0.5)) +
      
      scale_fill_brewer(palette = "Spectral")
      
    print(hb)
    
    if(printPDF) garbage <- dev.off()
    if(verbose) message(" done ")
  }
  
  # peptide detection (using modified.sequence)
  if (plotPEPDETECT) {
    if(verbose) message("--- Plot Peptide Detection (using modified.sequence)", appendLF = FALSE)
    if(printPDF) pdf(
      paste0(output_name,'.qcplot.PeptideDetection.pdf'),
      width = 10, #nsamples * 3
      height = 6,
      onefile = TRUE
    )
    peptDX <- peptCV[, .N, by = list(condition, pDX)]
    peptDXT <- peptDX[, list(N.total = sum(N)), by = list(condition)]
    peptDX <- merge(peptDX, peptDXT, by = "condition")
    peptDX$fxDx <- 100 * (peptDX$N / peptDX$N.total)
    
    ia <- ggplot(peptDX, aes(x = condition,
                             y = fxDx,
                             fill = factor(pDX))) +
      geom_bar(stat = "identity", alpha = 0.7, na.rm = TRUE) +
      ggrepel::geom_text_repel(
        aes(label = round(fxDx, digits = 1)),
        vjust = 1,
        size = 2,
        position = position_stack()
      ) +
      xlab("Condition") + ylab("Counts") +
      # ggtitle("Frequency of peptides detection") +
      labs(title = "Frequency of peptides detection",
           fill = "n") +
      theme_linedraw() +
      theme(legend.text = element_text(size = 8)) +
      theme(axis.text.x = element_text(angle = 90, 
                                       size = 8)) +
      theme(axis.text.y = element_text(size = 10)) +
      theme(plot.title = element_text(size = 12)) +
      scale_fill_brewer(palette = "Set1")
      
    print(ia)
    
    if(printPDF) garbage <- dev.off()
    if(verbose) message(" done ")
  }
  
  
  #m Protein Intensity CV
  if (plotPROTICV) {
    if(verbose) message("--- Plot Protein Intensity CV", appendLF = FALSE)
    if(printPDF) pdf(
      paste0(output_name,'.qcplot.ProteinIntensityCV.pdf'),
      width = 10, #nsamples * 3
      height = 6,
      onefile = TRUE
    )
    protCV <- data.table::data.table(subset(evidencekeys,!is.na(intensity)))
    protCV$intensity <- as.numeric(protCV$intensity)
    protCV <- protCV[, list(intensity = sum(intensity, na.rm = TRUE)), by = list(condition, bioreplicate, proteins)]
    protCV <- protCV[, list(pCV = 100 * (sd(intensity) / mean(intensity)),
                            sumInt = sum(intensity),
                            pDX = length(unique(bioreplicate))), 
                     by = list(condition, proteins)]
      
      
    protCV <- protCV[, bin.all := .artms_qcut(sumInt, 4)]
    protCV <- protCV[, bin.condition := .artms_qcut(sumInt, 4), by = condition]
    
    ja <- ggplot(subset(protCV,!is.na(pCV)), aes(condition, pCV)) +
      geom_boxplot(varwidth = TRUE, aes(fill = factor(condition)), alpha = 0.7, na.rm = TRUE) +
      geom_text(
        data = aggregate(pCV ~ condition, subset(protCV,!is.na(pCV)), median),
        aes(label = round(pCV, digits = 1),
            y = max(protCV$pCV, na.rm = TRUE) + 1),
        size = 2) +
      xlab("Condition") + ylab("Coefficient of variance (%)") +
      labs(title = "Distribution of Protein intensity CV",
           fill = "") +
      theme_linedraw() +
      theme(legend.text = element_text(size = 8)) +
      theme(axis.text.x = element_text(angle = 90, 
                                       size = 8)) +
      theme(axis.text.y = element_text(size = 10)) +
      theme(axis.title.y = element_text(size = 10)) +
      theme(plot.title = element_text(size = 12)) +
      scale_fill_brewer(palette = "Spectral")
      
    print(ja)
    
    jb <- ggplot(subset(protCV,!is.na(pCV)), aes(interaction(condition, bin.condition), pCV)) +
      geom_boxplot(varwidth = TRUE, aes(fill = factor(condition)), alpha = 0.7, na.rm = TRUE) +
      geom_text(
        data = aggregate(
          pCV ~ condition + bin.condition,
          subset(protCV,!is.na(pCV)),
          median
        ),
        aes(
          label = round(pCV, digits = 1),
          y = max(protCV$pCV, na.rm = TRUE) + 1
        ),
        size = 2,
        angle = 0
      ) +
      geom_text(
        data = aggregate(
          pCV ~ condition + bin.condition,
          subset(protCV,!is.na(pCV)),
          length
        ),
        aes(label = round(pCV, digits = 1), y = 0),
        size = 2,
        angle = 0
      ) +
      xlab("Condition") + ylab("Coefficient of variance (%)") +
      labs(title = "Distribution of Protein intensity CV",
           caption = "Proteins ranked by summed intensity and CV calculated for each protein.
           Each condition shows 4 distribution (box) for low (1) to high (4) intensity proteins.
           Top: overall median CV within each condition
           Bottom: number of protein groups used to calculate CV",
           fill = "") +
      theme_linedraw() +
      theme(legend.text = element_text(size = 8)) +
      theme(axis.text.x = element_text(angle = 90, 
                                       size = 8)) +
      theme(axis.text.y = element_text(size = 10)) +
      theme(axis.title.y = element_text(size = 10)) +
      theme(plot.title = element_text(size = 12)) +
      scale_fill_brewer(palette = "Spectral")
      
    print(jb)
    
    if(printPDF) garbage <- dev.off()
    if(verbose) message(" done ")
  }
  
  # Protein detection
  if (plotPROTDETECT) {
    if(verbose) message("--- Plot Protein Detection", appendLF = FALSE)
    if(printPDF) pdf(
      paste0(output_name,'.qcplot.ProteinDetection.pdf'),
      width = 10, #nsamples * 3
      height = 6,
      onefile = TRUE
    )
    protDX <- protCV[, .N, by = list(condition, pDX)]
    protDXT <- protDX[, list(N.total = sum(N)), by = list(condition)]
    protDX <- merge(protDX, protDXT, by = "condition")
    protDX$fxDx <- 100 * (protDX$N / protDX$N.total)
    
    ka <- ggplot(protDX, aes(
      x = condition,
      y = fxDx,
      fill = factor(pDX)
    )) +
      geom_bar(stat = "identity", alpha = 0.7, na.rm = TRUE) +
      ggrepel::geom_text_repel(
        aes(label = round(fxDx, digits = 1)),
        vjust = 1,
        size = 2,
        position = position_stack()
      ) +
      xlab("Condition") + ylab("Counts") +
      labs(title = "Protein Detection Frequency across bioreplicates by condition",
           fill = "n") +
      theme_linedraw() +
      theme(legend.text = element_text(size = 8)) +
      theme(axis.text.x = element_text(angle = 90, 
                                       size = 8)) +
      theme(axis.text.y = element_text(size = 10)) +
      scale_fill_brewer(palette = "Set1")
      
    print(ka)
    
    kb <- ggplot(subset(evidencekeys,!is.na(intensity)), aes(bioreplicate, log2(intensity))) +
      geom_boxplot(varwidth = TRUE,
                   aes(fill = potential.contaminant),
                   alpha = 0.7, na.rm = TRUE) +
      #ggrepel::geom_text_repel(data = aggregate(intensity ~ bioreplicate + potential.contaminant, subset(evidencekeys, !is.na(intensity)), median), aes(label = round(log2(intensity), digits=1), y = log2(max(evidencekeys$intensity, na.rm= TRUE))+0.5 ), size = 15) +
      xlab("Experiment") + 
      ylab("Log2 Intensity") +
      labs(title = "Peptide feature intensity distribution",
           fill = "Contaminants") +
      theme_linedraw() +
      theme(legend.text = element_text(size = 8)) +
      theme(axis.text.x = element_text(angle = 90, 
                                       size = 8)) +
      theme(axis.text.y = element_text(size = 10)) +
      theme(axis.title.y = element_text(size = 10)) +
      theme(plot.title = element_text(size = 12)) +
      scale_fill_brewer(palette = "Set1")
      
    print(kb)
    
    if (isFRACTION) {
      kc <-
        ggplot(subset(evidencekeys,!is.na(intensity)), aes(bioreplicate, log2(intensity))) +
        geom_boxplot(varwidth = TRUE,
                     aes(fill = potential.contaminant),
                     alpha = 0.7, na.rm = TRUE) +
        #ggrepel::geom_text_repel(data = aggregate(intensity ~ bioreplicate + potential.contaminant, subset(evidencekeys, !is.na(intensity)), median), aes(label = round(log2(intensity), digits=1), y = log2(max(evidencekeys$intensity, na.rm= TRUE))+0.5 ), size = 15) +
        xlab("Experiment") + ylab("Log2 Intensity") +
        facet_wrap( ~ fraction, ncol = 5) +
        ggtitle("Peptide feature intensity distribution by Fraction") +
        theme(legend.text = element_text(size = 8)) +
        theme(axis.text.x = element_text(angle = 90, 
                                         size = 10)) +
        theme(axis.text.y = element_text(size = 10)) +
        theme(axis.title.y = element_text(size = 10)) +
        theme(plot.title = element_text(size = 12)) +
        scale_fill_brewer(palette = "Set1")
      print(kc)
      
    }
    
    kd <- ggplot(evidencekeys, aes(x = log2(intensity))) +
      geom_density(alpha = .5, aes(fill = type), na.rm = TRUE) +
      labs(title = "Peptide feature intensity distribution by ID Type") +
      facet_wrap( ~ bioreplicate, ncol = 5) +
      theme_linedraw() +
      theme(legend.text = element_text(size = 8)) +
      theme(axis.text.x = element_text(size = 8)) +
      theme(axis.text.y = element_text(size = 10)) +
      theme(axis.title.y = element_text(size = 10)) +
      theme(plot.title = element_text(size = 12)) +
      scale_fill_brewer(palette = "Set1")
    
    print(kd)
    
    if(printPDF) garbage <- dev.off()
    if(verbose) message(" done ")
  }
  
  
  if (plotIDoverlap) {
    if(verbose) message("--- Plot ID overlap", appendLF = FALSE)
    if(printPDF) pdf(
      paste0(output_name,'.qcplot.ID-Overlap.pdf'),
      width = 20, #nsamples * 3
      height = 20,
      onefile = TRUE
    )
    ovlSeq <- subset(
      evidencekeys,!is.na(intensity) & !is.na(ms.ms.count),
      select = c("bioreplicate", "sequence")
    )
    ovlSeq <- unique(ovlSeq)
    ovlSeq <- .artms_sampleOverlap(ovlSeq,
                                   sampleID = 'bioreplicate',
                                   referenceID = 'sequence')
    ovlSeqM <- ovlSeq$M
    ovlSeqM[ovlSeqM == 1] <- NA
    # coordinates to move the keys
    lmat = rbind(c(0,3),c(2,1),c(0,4))
    lwid = c(1.5,4)
    lhei = c(1.5,4,1)
    gplots::heatmap.2(
      ovlSeqM * 100,
      col = hmcol,
      Rowv = FALSE,
      Colv = "Rowv",
      main = "Pairwise peptide identification overlap",
      cexRow = 1,
      cexCol = 1,
      margins = c(20, 20),
      cellnote = round(ovlSeqM * 100, 1),
      notecex = 1,
      notecol = "black",
      trace = "none",
      key = FALSE,
      keysize = 1,
      scale = "none",
      symm = TRUE,
      colsep = c(seq_len(length(
        colnames(ovlSeqM)
      ))),
      rowsep = c(seq_len(length(
        rownames(ovlSeqM)
      ))),
      sepcolor = "white",
      sepwidth = c(0.01, 0.01),
      dendrogram = "none",
      density.info = "none"
    )
    
    if (isFRACTION) {
      ovlSeqFx <- subset(
        evidencekeys,!is.na(intensity) & !is.na(ms.ms.count),
        select = c("bioreplicate", "fraction", "sequence")
      )
      ovlSeqFx <- unique(ovlSeqFx)
      ovlSeqFx$bioreplicate <- paste(ovlSeqFx$fraction,
                                     ovlSeqFx$bioreplicate,
                                     sep = ".")
      ovlSeqFx <- .artms_sampleOverlap(ovlSeqFx,
                                       sampleID = 'bioreplicate',
                                       referenceID = 'sequence')
      ovlSeqFxM <- ovlSeqFx$M
      ovlSeqFxM[ovlSeqFxM == 1] <- NA
      gplots::heatmap.2(
        ovlSeqFxM * 100,
        col = hmcol,
        Rowv = FALSE,
        Colv = "Rowv",
        cexRow = 0.5,
        cexCol = 0.5,
        margins = c(40, 40),
        #cellnote=round(ovlSeqFxM*100,1),
        #notecex=3,
        main = "Pairwise peptide identification overlap by Fraction",
        #notecol="black",
        trace = "none",
        key = TRUE,
        keysize = 1,
        scale = "none",
        symm = TRUE,
        colsep = c(seq_len(length(
          colnames(ovlSeqFxM)
        ))),
        rowsep = c(seq_len(length(
          rownames(ovlSeqFxM)
        ))),
        sepcolor = "white",
        sepwidth = c(0.0001, 0.0001),
        dendrogram = "none"
      )
    }
    
    
    ovlProt <- subset(
      evidencekeys,!is.na(intensity) & !is.na(ms.ms.count),
      select = c("bioreplicate", "proteins")
    )
    ovlProt <- unique(ovlProt)
    ovlProt <- .artms_sampleOverlap(ovlProt,
                                    sampleID = 'bioreplicate',
                                    referenceID = 'proteins')
    ovlProtM <- ovlProt$M
    ovlProtM[ovlProtM == 1] <- NA
    gplots::heatmap.2(
      ovlProtM * 100,
      col = hmcol,
      Rowv = FALSE,
      Colv = "Rowv",
      cexRow = 1,
      cexCol = 1,
      margins = c(40, 40),
      cellnote = round(ovlProtM * 100, 1),
      notecex = 1,
      main = "Pairwise protein identification overlap",
      notecol = "black",
      trace = "none",
      key = FALSE,
      keysize = 1,
      scale = "none",
      symm = TRUE,
      colsep = c(seq_len(length(
        colnames(ovlProtM)
      ))),
      rowsep = c(seq_len(length(
        rownames(ovlProtM)
      ))),
      sepcolor = "white",
      sepwidth = c(0.01, 0.01),
      dendrogram = "none"
    )
    if(printPDF) garbage <- dev.off()
    if(verbose) message(" done ")
  }
  
  
  if (plotPCA) {
    if(verbose) message("--- Plot PCA and Inter-Correlation (WARNING: it might take a long time. Please, be patient)")
    if(printPDF) pdf(
      paste0(output_name,'.qcplot.PCA.pdf'),
      width = 12, #nsamples * 3
      height = 10,
      onefile = TRUE
    )
    
    peptintmtx <- subset(evidencekeys,
                         !is.na(intensity),
                         select = c("bioreplicate", "sequence", "intensity"))

    peptintmtx$intensity <- as.integer64(peptintmtx$intensity)
  
    ##LEGACY
    # peptintmtx <- data.table::dcast(peptintmtx, 
    #                                 sequence ~ bioreplicate, 
    #                                 sum, 
    #                                 value.var = "intensity")
    
    peptintmtx <- peptintmtx %>% 
      tidyr::pivot_wider(id_cols = sequence, 
                         names_from = bioreplicate, 
                         values_from = intensity, 
                         values_fn = list(intensity = sum))
    
    
    peptintmtx <- as.matrix(peptintmtx[, -1])
    peptintmtx <- log2(peptintmtx)
    peptintmtx[!is.finite(peptintmtx)] <- NA
    
    if(nsamples < 20){
      pairs(peptintmtx,
            upper.panel = .artms_panelCor,
            diag.panel = .artms_panelHist,
            lower.panel = panel.smooth,
            pch = ".")
    }else{
      message("\t(-) Skip peptide-based correlation matrix (too many samples)")
    }
    
    mpept <- peptintmtx
    mpept[is.na(mpept)] <- 0
    pept.pca <- prcomp(t(mpept))
    pept.pcax <- as.data.frame(pept.pca$x)
    pept.pcax <- merge(unique(evidencekeys[, c("bioreplicate", "condition")]),
            pept.pcax,
            by.x = "bioreplicate",
            by.y = "row.names")
    
    ma <- ggplot(pept.pcax, aes(PC1, PC2, fill = condition)) + #, color=condition, shape=condition
      geom_point(alpha = .8, size = 5, na.rm = TRUE) +
      geom_label_repel(
        aes(label = bioreplicate),
        # label.size = 1,
        box.padding   = 0.8,
        point.padding = 0.2
        # segment.color = 'grey50'
      ) +
      labs(title = "Principal Component Analysis (Peptide Intensity)",
           fill = "") +
      theme_linedraw()

    print(ma)
    
    protintmtx <- subset(
      evidencekeys,
      !is.na(intensity),
      select = c("bioreplicate", "proteins", "intensity")
    )
    
    protintmtx$intensity <- as.integer64(protintmtx$intensity)
    ##LEGACY
    # protintmtx <- data.table::dcast(protintmtx, 
    #                                 proteins ~ bioreplicate, 
    #                                 sum, 
    #                                 value.var = "intensity")

    protintmtx <- protintmtx %>% 
      tidyr::pivot_wider(id_cols = proteins, 
                         names_from = bioreplicate, 
                         values_from = intensity, 
                         values_fn = list(intensity = sum ))
    
    
    protintmtx <- as.matrix(protintmtx[, -1])
    protintmtx <- log2(protintmtx)
    protintmtx[!is.finite(protintmtx)] <- NA
    
    if(nsamples < 20){
      pairs(
        protintmtx,
        upper.panel = .artms_panelCor,
        diag.panel = .artms_panelHist,
        lower.panel = panel.smooth,
        pch = "."
      )
    }else{
      message("\t(-) Skip Protein-based correlation matrix (too many samples)")
    }

    mprot <- protintmtx
    mprot[is.na(mprot)] <- 0
    prot.pca <- prcomp(t(mprot))
    prot.pcax <- as.data.frame(prot.pca$x)
    prot.pcax <- merge(unique(evidencekeys[, c("bioreplicate", "condition")]),
                       prot.pcax,
                       by.x = "bioreplicate",
                       by.y = "row.names")
      
    
    mb <- ggplot(prot.pcax, aes(PC1, PC2, fill = condition)) + #, color=condition, shape=condition
      geom_point(alpha = .8, size = 5, na.rm = TRUE) +
      geom_label_repel(
        aes(label = bioreplicate),
        #label.size = 2,
        box.padding   = 0.8,
        point.padding = 0.1,
        segment.color = 'grey50'
      ) +
      labs(title = "Principal Component Analysis (Protein Intensity)",
           fill = "") +
      theme_linedraw()
    
    print(mb)
    
    if(printPDF) garbage <- dev.off()
  }
  
  
  if (plotSP) {

    if(verbose) message("--- Plot Sample Preparation", appendLF = FALSE)
    if(printPDF) pdf(
      paste0(output_name,'.qcplot.SamplePrep.pdf'),
      width = 10, #nsamples * 3
      height = 6,
      onefile = TRUE
    )
    
    if("missed.cleavages" %in% colnames(evidencekeys)){
      evidence.misscleavages <- data.table(subset(evidencekeys, ms.ms.count > 0))
      misscleavages.tot <- evidence.misscleavages[, .N, by = list(bioreplicate)]
      misscleavages.dt <- evidence.misscleavages[, .N,
                                                 by = list(bioreplicate, missed.cleavages)]
      misscleavages.dt <- merge(misscleavages.dt,
                                misscleavages.tot,
                                by = "bioreplicate",
                                all = TRUE)
      misscleavages.dt$mc <- as.numeric(format(100 * (misscleavages.dt$N.x / misscleavages.dt$N.y),
                                               digits = 5))
      naplot <- ggplot(misscleavages.dt, aes(x = bioreplicate,
                                         y = mc,
                                         fill = as.factor(missed.cleavages),
                                         label = mc)) +
        geom_bar(stat = "identity", alpha = 0.7, na.rm = TRUE) +
        geom_text(
          aes(label = round(mc, digits = 1)),
          vjust = 0,
          size = 2,
          position = position_stack()
        ) +
        xlab("Experiment") + ylab("Fraction (percentage)") +
        labs(title = "Missing cleavages",
             fill = "") +
        theme_linedraw() +
        theme(legend.text = element_text(size = 8)) +
        theme(axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          size = 8
        )) +
        theme(axis.text.y = element_text(size = 10)) +
        theme(axis.title.y = element_text(size = 10)) +
        theme(plot.title = element_text(size = 12)) +
        scale_fill_brewer(palette = "Set1")
        
      print(naplot)
    }


    oxMet.df <-
      plyr::ddply(
        evidencekeys,
        c("bioreplicate", "condition"),
        plyr::summarise,
        pct.OxM = (length(oxidation..m.[oxidation..m. > 0]) / length(sequence)) * 100)
    
    nb <- ggplot(oxMet.df,
                 aes(
                   x = bioreplicate,
                   y = pct.OxM,
                   fill = condition,
                   label = pct.OxM
                 )) +
      geom_bar(stat = "identity", alpha = 0.7, na.rm = TRUE) +
      geom_text(
        aes(label = round(pct.OxM,
                          digits = 1)),
        hjust = 1,
        size = 2.7,
        angle = 0,
        position = position_dodge(width = 1)
      ) +
      xlab("Experiment") +
      ylab("Fraction (percentage)") +
      ggtitle("Percentage of peptides with at least 1 Methionine oxidized") +
      theme_linedraw() +
      theme(legend.text = element_text(size = 8)) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 8
      )) +
      theme(axis.text.y = element_text(size = 10)) +
      theme(axis.title.y = element_text(size = 10)) +
      theme(plot.title = element_text(size = 12)) +
      scale_fill_brewer(palette = "Spectral")
    
    print(nb)
    
    if(printPDF) garbage <- dev.off()
    if(verbose) message("... done")
  }
  
  if(verbose) message(">> QC extended completed")

} # END OF artmsQualityControlEvidenceExtended  


# @title Quantiles
#
# @description Quantiles
# @param x (datatype) Describe
# @param n (datatype) Describe
# @return What does it return?
# @keywords internal, stats
.artms_qcut <- function(x, n) {
  quantiles <- seq(0, 1, length.out = n + 1)
  cutpoints <- unname(quantile(x, quantiles, na.rm = TRUE))
  as.character(cut(
    x,
    cutpoints,
    labels = seq_len(n),
    include.lowest = TRUE
  ))
}

# @title Sample overlap
#
# @description Describe
# @param data (datatype) Describe
# @param sampleID (datatype) Describe
# @param referenceID (datatype) Describe
# @return What does it return?
# @keywords internal, stats
.artms_sampleOverlap <- function(data,
                                 sampleID = 'bioreplicate',
                                 referenceID = 'Sequence') {
  pept <- by(data,
             data[, sampleID],
             function(x)
               unique(as.character(x[, referenceID])))
  pepmatch <- lapply(pept,
                     function(x, y) {
                       lapply(y, function(y)
                         length(intersect(x, y)) / length(union(x, y)))
                     }, pept)
  m <- matrix(unlist(pepmatch),
              nrow = length(pept),
              ncol = length(pept),
              byrow = TRUE)
    
  colnames(m) <- rownames(m) <- names(pept)
  ans <- list(Peptides = pept, M = m)
  return(ans)
}

# @title Panel Cor
#
# @description Describe
# @param x (datatype) Describe
# @param y (datatype) Describe
# @param digits A number representing...
# @param prefix A prefix for the...
# @param cex.cor A Cex cor....
# @param ... No idea what this does
# @return What does it return?
# @keywords internal, stats
.artms_panelCor <- function(x,
                            y,
                            digits = 2,
                            prefix = "",
                            cex.cor,
                            ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use = "pairwise.complete.obs", method = "pearson"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste(prefix, txt, sep = "")
  if (missing(cex.cor))
    cex.cor <- 0.8 / strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
  

# @title Panel Histogram
# @description Describe
# @param x (datatype) Describe
# @param ... No idea what this does
# @return What does it return?
# @keywords internal, stats
.artms_panelHist <- function(x, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(usr[seq_len(2)], 0, 1.5))
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks
  nB <- length(breaks)
  y <- h$counts
  y <- y / max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}