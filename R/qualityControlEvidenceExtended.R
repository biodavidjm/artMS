# ------------------------------------------------------------------------------
#' @title Extended Quality Control of the MaxQuant evidence.txt file
#'
#' @description Performs quality control based on the information available in
#' the MaxQuant `evidence.txt` file.
#' @param evidence_file (char or data.frame) The evidence file path and name, or
#' data.frame
#' @param keys_file (char or data.frame) The keys file path and name or
#' data.frame
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
#' @param plotPROTEINS (logical) `TRUE` generates protein statistics plot: 
#' Page 1 shows the number of protein groups confidently identified in each 
#' BioReplicate. If replicates are present, Page 2 shows the mean number of 
#' protein groups per condition with error bar showing the standard error of 
#' the mean. Note that potential contaminant proteins are plotted separately.
#' Pages 3 and 4 show peptide identification intersection between 
#' BioReplicates (the bars are ordered by degree or frequency, respectively), 
#' and Page 4 shows the intersections across conditions instead of 
#' BioReplicates.
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
#' @param plotIC (logical) `TRUE` generates pairwise intensity correlation: 
#' Page 1 and 3 show pairwise peptide and protein intensity correlation and 
#' scatter plot between any 2 BioReplicates, respectively. Page 2 and 4 show 
#' principal component analysis at the intensity level for both peptide and 
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
#' test <- artmsQualityControlEvidenceExtended  (evidence_file = NULL,
#' keys_file = NULL)
#' @export
artmsQualityControlEvidenceExtended   <- function(evidence_file,
                                                 keys_file,
                                                 plotPSM = TRUE,
                                                 plotIONS = TRUE,
                                                 plotTYPE = TRUE,
                                                 plotPEPTIDES = TRUE,
                                                 plotPROTEINS = TRUE,
                                                 plotPIO = TRUE,
                                                 plotCS = TRUE,
                                                 plotME = TRUE,
                                                 plotMOCD = TRUE,
                                                 plotPEPICV = TRUE,
                                                 plotPEPDETECT = TRUE,
                                                 plotPROTICV = TRUE,
                                                 plotPROTDETECT = TRUE,
                                                 plotIDoverlap = TRUE,
                                                 plotIC = TRUE,
                                                 plotSP = TRUE,
                                                 printPDF = TRUE,
                                                 verbose = TRUE) {
  if(verbose)
    message(">>EXTENDED QUALITY CONTROL ANALYSIS
(evidence.txt based)---------------------- ")
  
  if (is.null(evidence_file) & is.null(keys_file)) {
    return("You need to provide both evidence and keys")
  }
  
  if(any(missing(evidence_file) | 
         missing(keys_file)))
    stop("Missed (one or many) required argument(s)
         Please, check the help of this function to find out more")
  
  hmcol <- colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(256)
  
  # ----------------------------------------------------------------------------
  # DATA PREPARATION
  evidencekeys <-
    artmsMergeEvidenceAndKeys(evidence_file, 
                               keys_file,
                               verbose = verbose)
  colnames(evidencekeys) <- tolower(colnames(evidencekeys))
  
  if (any(grepl("+", evidencekeys$reverse))) {
    evidencekeys <- subset(evidencekeys, reverse != "+")
  }
  
  # Add potential.contaminant column in case search was performed without it
  `%ni%` <- Negate(`%in%`)
  if ("potential.contaminant" %ni% colnames(evidencekeys)) {
    evidencekeys$potential.contaminant <- ""
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
  
  if ("fraction" %in% colnames(evidencekeys)) {
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
  
  evidence3 <-
    plyr::ddply(
      evidence2,
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
  oversampling0 <-
    evidencekeys.dt[, .N, by = list(ms.ms.count, bioreplicate)]
  oversampling0$MSMS.counts <-
    findInterval(oversampling0$ms.ms.count, seq(1, 3, by = 1))
  oversampling0 <-
    oversampling0[, list(N = sum(N)), by = list(MSMS.counts, bioreplicate)]
  oversampling0$MSMS.counts <-
    paste("n=", oversampling0$MSMS.counts, sep = "")
  oversampling0$MSMS.counts <-
    gsub("n=3", "n=3+", oversampling0$MSMS.counts)
  oversampling0.total <-
    evidencekeys.dt[, .N, by = list(bioreplicate)]
  oversampling0 <-
    merge(oversampling0,
          oversampling0.total,
          by = "bioreplicate",
          all = TRUE)
  oversampling0$FxOverSamp <-
    as.numeric(format(100 * (oversampling0$N.x / oversampling0$N.y), digits = 3))
  
  oversampling1 <-
    evidencekeys.dt[intensity > 0, .N, by = list(ms.ms.count, bioreplicate)]
  oversampling1$MSMS.counts <-
    findInterval(oversampling1$ms.ms.count, seq(1, 3, by = 1))
  oversampling1 <-
    oversampling1[, list(N = sum(N)), by = list(MSMS.counts, bioreplicate)]
  oversampling1$MSMS.counts <-
    paste("n=", oversampling1$MSMS.counts, sep = "")
  oversampling1$MSMS.counts <-
    gsub("n=3", "n=3+", oversampling1$MSMS.counts)
  oversampling1.total <-
    evidencekeys.dt[intensity > 0, .N, by = list(bioreplicate)]
  oversampling1 <-
    merge(oversampling1,
          oversampling1.total,
          by = "bioreplicate",
          all = TRUE)
  oversampling1$FxOverSamp <-
    as.numeric(format(100 * (oversampling1$N.x / oversampling1$N.y), digits = 3))
  
  oversampling2 <-
    evidencekeys.dt[intensity > 0 &
                      ms.ms.count > 0, .N, by = list(ms.ms.count, bioreplicate)]
  oversampling2$MSMS.counts <-
    findInterval(oversampling2$ms.ms.count, seq(1, 3, by = 1))
  oversampling2 <-
    oversampling2[, list(N = sum(N)), by = list(MSMS.counts, bioreplicate)]
  oversampling2$MSMS.counts <-
    paste("n=", oversampling2$MSMS.counts, sep = "")
  oversampling2$MSMS.counts <-
    gsub("n=3", "n=3+", oversampling2$MSMS.counts)
  oversampling2.total <-
    evidencekeys.dt[intensity > 0 &
                      ms.ms.count > 0, .N, by = list(bioreplicate)]
  oversampling2 <-
    merge(oversampling2,
          oversampling2.total,
          by = "bioreplicate",
          all = TRUE)
  oversampling2$FxOverSamp <-
    as.numeric(format(100 * (oversampling2$N.x / oversampling2$N.y), digits = 3))
  
  
  # Charge state
  chargeState <-
    evidencekeys.dt[, .N, by = list(charge, bioreplicate)]
  chargeState$charge <- paste("z=", chargeState$charge, sep = "")
  chargeState.total <- evidencekeys.dt[, .N, by = list(bioreplicate)]
  chargeState <-
    merge(chargeState,
          chargeState.total,
          by = "bioreplicate",
          all = TRUE)
  chargeState$FxOverSamp <-
    as.numeric(format(100 * (chargeState$N.x / chargeState$N.y), digits = 1))
  
  chargeStateCond <-
    evidencekeys.dt[, .N, by = list(charge, condition)]
  chargeStateCond$charge <-
    paste("z=", chargeStateCond$charge, sep = "")
  chargeStateCond.total <- evidencekeys.dt[, .N, by = list(condition)]
  chargeStateCond <-
    merge(chargeStateCond,
          chargeStateCond.total,
          by = "condition",
          all = TRUE)
  chargeStateCond$FxOverSamp <-
    as.numeric(format(100 * (
      chargeStateCond$N.x / chargeStateCond$N.y
    ), digits = 1))
  
  # Type
  mstype <-
    evidencekeys.dt[, .N, by = list(type, bioreplicate, condition)]
  mstype.total <- evidencekeys.dt[, .N, by = list(bioreplicate)]
  mstype <-
    merge(mstype, mstype.total, by = "bioreplicate", all = TRUE)
  mstype$FxOverSamp <-
    as.numeric(format(100 * (mstype$N.x / mstype$N.y), digits = 2))
  
  # ----------------------------------------------------------------------------
  # PLOTS
  if(verbose) message(">> GENERATING QC PLOTS ")
  
  nsamples <- length(unique(evidencekeys$bioreplicate))
  nconditions <- length(unique(evidencekeys$condition))
  
  # Check the largest number so that width and height of pdf is not too large
  if (nsamples > 20) {
    nsamples <- 20
  }
  
  if (nconditions > 7) {
    nconditions <- 7
  }
  
  ### PSM
  if (plotPSM) {
    if(verbose) message("--- Plot PSM")
    if(printPDF) pdf(
      'QC_Plots_PSM.pdf',
      width = nsamples * 3,
      height = 20,
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
               alpha = 0.7) +
      geom_text(
        aes(label = round(PSMs, digits = 0)),
        hjust = 0.5,
        vjust = -0.5,
        size = 10,
        position = position_dodge(width = 1)
      ) +
      facet_wrap( ~ potential.contaminant, ncol = 1) +
      xlab("Experiment") + ylab("Counts") +
      ggtitle("Number of PSMs: bottom = Potential contaminants; top = non-contaminants") +
      theme(legend.text = element_text(size = 20)) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 0,
        size = 20
      )) +
      theme(axis.text.y = element_text(size = 20)) +
      theme(axis.title.x = element_text(size = 30)) +
      theme(axis.title.y = element_text(size = 30)) +
      theme(plot.title = element_text(size = 40)) +
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
               alpha = 0.7) +
      geom_errorbar(
        aes(ymin = PSMs.mean - PSMs.sem, ymax = PSMs.mean + PSMs.sem),
        width = .2,
        position = position_dodge(.9)
      ) +
      geom_text(
        aes(label = round(PSMs.mean, digits = 0)),
        hjust = 0.5,
        vjust = -1.5,
        size = 10,
        position = position_dodge(width = 1)
      ) +
      xlab("Condition") + ylab("Counts") +
      ggtitle(
        "Mean number of PSMs per condition for contaminants and non-contaminants, error bar= std error of the mean"
      ) +
      theme(legend.text = element_text(size = 20)) +
      theme(axis.text.x = element_text(
        angle = 0,
        hjust = 1,
        size = 20
      )) +
      theme(axis.text.y = element_text(size = 20)) +
      theme(axis.title.x = element_text(size = 30)) +
      theme(axis.title.y = element_text(size = 30)) +
      theme(plot.title = element_text(size = 40)) +
      scale_fill_brewer(palette = "Spectral")
    print(ab)
    
    if ("fraction" %in% colnames(evidencekeys)) {
      ac <-
        ggplot(evidence2fx, aes(
          x = bioreplicate,
          y = PSMs,
          fill = factor(fraction)
        )) +
        geom_bar(stat = "identity", alpha = 0.7) +
        geom_text(
          aes(label = round(PSMs, digits = 0)),
          hjust = 0.5,
          vjust = 1.5,
          size = 7,
          position = position_stack()
        ) +
        facet_wrap( ~ potential.contaminant,
                    ncol = 1,
                    scales = "free") +
        xlab("Experiment") + ylab("Counts") +
        ggtitle(
          "Number of PSMs per Fraction: bottom = Potential contaminants; top = non-contaminants"
        ) +
        theme(legend.text = element_text(size = 20)) +
        theme(axis.text.x = element_text(
          angle = 90,
          hjust = 0,
          size = 20
        )) +
        theme(axis.text.y = element_text(size = 20)) +
        theme(axis.title.x = element_text(size = 30)) +
        theme(axis.title.y = element_text(size = 30)) +
        theme(plot.title = element_text(size = 40))
      print(ac)
    }
    
    if(printPDF) garbage <- dev.off()
    if(verbose) message(" done ")
  }
  
  
  ### IONS
  if (plotIONS) {
    if(verbose) message("--- Plot IONS")
    
    if(printPDF) pdf(
      'QC_Plots_IONS.pdf',
      width = nsamples * 3,
      height = 20,
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
               alpha = 0.7) +
      geom_text(
        aes(label = round(Ions, digits = 0)),
        hjust = 0.5,
        vjust = -0.5,
        size = 10,
        position = position_dodge(width = 1)
      ) +
      facet_wrap( ~ potential.contaminant, ncol = 1) +
      xlab("Experiment") + ylab("Counts") +
      ggtitle(
        "Number of unique Peptide Ions: bottom = Potential contaminants; top = non-contaminants"
      ) +
      theme(legend.text = element_text(size = 20)) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 20
      )) +
      theme(axis.text.y = element_text(size = 20)) +
      theme(axis.title.x = element_text(size = 30)) +
      theme(axis.title.y = element_text(size = 30)) +
      theme(plot.title = element_text(size = 40)) +
      scale_fill_brewer(palette = "Spectral")
    print(ba)
    
    bb <-
      ggplot(evidence3, aes(
        x = condition,
        y = Ions.mean,
        fill = factor(potential.contaminant)
      )) +
      geom_bar(stat = "identity",
               position = position_dodge(width = 1),
               alpha = 0.7) +
      geom_errorbar(
        aes(ymin = Ions.mean - Ions.sem, ymax = Ions.mean + Ions.sem),
        width = .2,
        position = position_dodge(.9)
      ) +
      geom_text(
        aes(label = round(Ions.mean, digits = 0)),
        hjust = 0.5,
        vjust = -0.5,
        size = 10,
        position = position_dodge(width = 1)
      ) +
      xlab("Condition") + ylab("Counts") +
      ggtitle(
        "Mean number of unique Peptide Ions per condition for contaminants (blue) and non-contaminants (red), error bar= std error of the mean"
      ) +
      theme(legend.text = element_text(size = 20)) +
      theme(axis.text.x = element_text(
        angle = 0,
        hjust = 1,
        size = 20
      )) +
      theme(axis.text.y = element_text(size = 20)) +
      theme(axis.title.x = element_text(size = 30)) +
      theme(axis.title.y = element_text(size = 30)) +
      theme(plot.title = element_text(size = 40)) +
      scale_fill_brewer(palette = "Set1")
    print(bb)
    
    if ("fraction" %in% colnames(evidencekeys)) {
      bc <-
        ggplot(evidence2fx, aes(
          x = bioreplicate,
          y = Ions,
          fill = factor(fraction)
        )) +
        geom_bar(stat = "identity", alpha = 0.7) +
        geom_text(
          aes(label = round(Ions, digits = 0)),
          hjust = 0.5,
          vjust = 1.5,
          size = 7,
          position = position_stack()
        ) +
        facet_wrap( ~ potential.contaminant,
                    ncol = 1,
                    scales = "free") +
        xlab("Experiment") + ylab("Counts") +
        ggtitle(
          "Number of unique Peptide Ions in each Fraction: bottom = Potential contaminants; top = non-contaminants"
        ) +
        theme(legend.text = element_text(size = 20)) +
        theme(axis.text.x = element_text(
          angle = 90,
          hjust = 0,
          size = 20
        )) +
        theme(axis.text.y = element_text(size = 20)) +
        theme(axis.title.x = element_text(size = 30)) +
        theme(axis.title.y = element_text(size = 30)) +
        theme(plot.title = element_text(size = 40))
      print(bc)
    }
    
    if(printPDF) garbage <- dev.off()
    
    if(verbose) message(" done ")
  }
  
  #### TYPE ######
  if (plotTYPE) {
    if(verbose) message("--- Plot TYPE")
    
    if(printPDF) pdf(
      'QC_Plots_TYPE.pdf',
      width = nsamples * 3,
      height = 20,
      onefile = TRUE
    )
    
    ca <-
      ggplot(mstype, aes(
        x = bioreplicate,
        y = FxOverSamp,
        fill = factor(type)
      )) +
      geom_bar(stat = "identity",
               position = position_stack(),
               alpha = 0.7) +
      ggrepel::geom_text_repel(
        aes(label = round(FxOverSamp, digits = 1)),
        vjust = 1,
        size = 10,
        position = position_stack()
      ) +
      xlab("Experiment") + ylab("Fraction") +
      ggtitle("Type of identification (MaxQuant type column)") +
      theme(legend.text = element_text(size = 20)) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 20
      )) +
      theme(axis.text.y = element_text(size = 20)) +
      theme(axis.title.x = element_text(size = 30)) +
      theme(axis.title.y = element_text(size = 30)) +
      theme(plot.title = element_text(size = 40)) +
      scale_fill_brewer(palette = "Set1")
    print(ca)
    
    if(printPDF) garbage <- dev.off()
    
    if(verbose) message(" done ")
  }
  
  
  
  ### PEPTIDES ###
  if (plotPEPTIDES) {
    if(verbose) message("--- Plot PEPTIDES")
    
    if(printPDF) pdf(
      'QC_Plots_PEPTIDES.pdf',
      width = nsamples * 3,
      height = 20,
      onefile = TRUE
    )
    
    da <-
      ggplot(evidence2, aes(
        x = bioreplicate,
        y = Peptides,
        fill = factor(condition)
      )) +
      geom_bar(stat = "identity",
               position = position_dodge(width = 1),
               alpha = 0.7) +
      geom_text(
        aes(label = round(Peptides, digits = 0)),
        hjust = 0.5,
        vjust = -0.5,
        size = 10,
        position = position_dodge(width = 1)
      ) +
      facet_wrap( ~ potential.contaminant, ncol = 1) +
      xlab("Experiment") + ylab("Counts") +
      ggtitle("Number of unique Peptides: bottom = Potential contaminants; top = non-contaminants") +
      theme(legend.text = element_text(size = 20)) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 20
      )) +
      theme(axis.text.y = element_text(size = 20)) +
      theme(axis.title.x = element_text(size = 30)) +
      theme(axis.title.y = element_text(size = 30)) +
      theme(plot.title = element_text(size = 40)) +
      scale_fill_brewer(palette = "Spectral")
    print(da)
    
    db <-
      ggplot(evidence3,
             aes(
               x = condition,
               y = Peptides.mean,
               fill = factor(potential.contaminant)
             )) +
      geom_bar(stat = "identity",
               position = position_dodge(width = 1),
               alpha = 0.7) +
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
        size = 10,
        position = position_dodge(width = 1)
      ) +
      xlab("Condition") + ylab("Counts") +
      ggtitle(
        "Mean number of unique Peptides per condition for contaminants (blue) and non-contaminants (red), error bar= std error of the mean"
      ) +
      theme(legend.text = element_text(size = 20)) +
      theme(axis.text.x = element_text(
        angle = 0,
        hjust = 1,
        size = 20
      )) +
      theme(axis.text.y = element_text(size = 20)) +
      theme(axis.title.x = element_text(size = 30)) +
      theme(axis.title.y = element_text(size = 30)) +
      theme(plot.title = element_text(size = 40)) +
      scale_fill_brewer(palette = "Set1")
    print(db)
    
    if ("fraction" %in% colnames(evidencekeys)) {
      dc <-
        ggplot(evidence2fx,
               aes(
                 x = bioreplicate,
                 y = Peptides,
                 fill = factor(fraction)
               )) +
        geom_bar(stat = "identity", alpha = 0.7) +
        geom_text(
          aes(label = round(Peptides, digits = 0)),
          hjust = 0.5,
          vjust = 1.5,
          size = 7,
          position = position_stack()
        ) +
        facet_wrap( ~ potential.contaminant,
                    ncol = 1,
                    scales = "free") +
        xlab("Experiment") + ylab("Counts") +
        ggtitle(
          "Number of unique Peptides in each Fraction: bottom = Potential contaminants; top = non-contaminants"
        ) +
        theme(legend.text = element_text(size = 20)) +
        theme(axis.text.x = element_text(
          angle = 90,
          hjust = 0,
          size = 20
        )) +
        theme(axis.text.y = element_text(size = 20)) +
        theme(axis.title.x = element_text(size = 30)) +
        theme(axis.title.y = element_text(size = 30)) +
        theme(plot.title = element_text(size = 40))
      print(dc)
    }
    
    lst.pepExp <- list()
    for (i in unique(evidencekeys$bioreplicate)) {
      lst.pepExp[[i]] <-
        unique(subset(evidencekeys, ms.ms.count > 0 &
                        bioreplicate == i)[, c("sequence")])
    }
    
    UpSetR::upset(
      fromList(lst.pepExp),
      nsets = length(unique(evidencekeys$bioreplicate)),
      nintersects = 50,
      order.by = "degree"
    )
    
    UpSetR::upset(
      fromList(lst.pepExp),
      nsets = length(unique(evidencekeys$bioreplicate)),
      nintersects = 50,
      order.by = "freq"
    )
    
    lst.pepCond <- list()
    for (i in unique(evidencekeys$condition)) {
      lst.pepCond[[i]] <- unique(subset(evidencekeys,
                                        ms.ms.count > 0 &
                                          condition == i)[, c("sequence")])
    }
    
    #UpSetR::upset(fromList(lst.pepCond), nsets=length(unique(evidencekeys$condition)), nintersects = 50, order.by = "degree")
    UpSetR::upset(
      fromList(lst.pepCond),
      nsets = length(unique(evidencekeys$condition)),
      nintersects = 50,
      order.by = "freq"
    )
    
    if(printPDF) garbage <- dev.off()
    if(verbose) message(" done ")
  }
  
  
  ## PROTEINS
  if (plotPROTEINS) {
    if(verbose) message("--- Plot PROTEINS")
    
    if(printPDF) pdf(
      'QC_Plots_PROTEINS.pdf',
      width = nsamples * 3,
      height = 20,
      onefile = TRUE
    )
    
    ea <-
      ggplot(evidence2, aes(
        x = bioreplicate,
        y = Proteins,
        fill = factor(condition)
      )) +
      geom_bar(stat = "identity",
               position = position_dodge(width = 1),
               alpha = 0.7) +
      geom_text(
        aes(label = round(Proteins, digits = 0)),
        hjust = 0.5,
        vjust = -0.5,
        size = 10,
        position = position_dodge(width = 1)
      ) +
      facet_wrap( ~ potential.contaminant, ncol = 1) +
      xlab("Experiment") + ylab("Counts") +
      ggtitle(
        "Number of unique Protein Groups: bottom = Potential contaminants; top = non-contaminants"
      ) +
      theme(legend.text = element_text(size = 20)) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 20
      )) +
      theme(axis.text.y = element_text(size = 20)) +
      theme(axis.title.x = element_text(size = 30)) +
      theme(axis.title.y = element_text(size = 30)) +
      theme(plot.title = element_text(size = 40)) +
      scale_fill_brewer(palette = "Spectral")
    print(ea)
    
    eb <-
      ggplot(evidence3,
             aes(
               x = condition,
               y = Proteins.mean,
               fill = factor(potential.contaminant)
             )) +
      geom_bar(stat = "identity",
               position = position_dodge(width = 1),
               alpha = 0.7) +
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
        size = 10,
        position = position_dodge(width = 1)
      ) +
      xlab("Condition") + ylab("Counts") +
      ggtitle(
        "Mean number of unique Proteins per condition for contaminants (blue) and non-contaminants (red), error bar= std error of the mean"
      ) +
      theme(legend.text = element_text(size = 20)) +
      theme(axis.text.x = element_text(
        angle = 0,
        hjust = 1,
        size = 20
      )) +
      theme(axis.text.y = element_text(size = 20)) +
      theme(axis.title.x = element_text(size = 30)) +
      theme(axis.title.y = element_text(size = 30)) +
      theme(plot.title = element_text(size = 40)) +
      scale_fill_brewer(palette = "Set1")
    print(eb)
    
    lst.protExp <- list()
    for (i in unique(evidencekeys$bioreplicate)) {
      lst.protExp[[i]] <-
        unique(subset(evidencekeys, ms.ms.count > 0 &
                        bioreplicate == i)[, c("proteins")])
    }
    
    UpSetR::upset(
      fromList(lst.protExp),
      nsets = length(unique(evidencekeys$bioreplicate)),
      nintersects = 50,
      order.by = "degree"
    )
    UpSetR::upset(
      fromList(lst.protExp),
      nsets = length(unique(evidencekeys$bioreplicate)),
      nintersects = 50,
      order.by = "freq"
    )
    
    lst.protCond <- list()
    for (i in unique(evidencekeys$condition)) {
      lst.protCond[[i]] <-
        unique(subset(evidencekeys, ms.ms.count > 0 &
                        condition == i)[, c("proteins")])
    }
    
    #UpSetR::upset(fromList(lst.protCond), nsets=length(unique(evidencekeys$condition)), nintersects = 50, order.by = "degree")
    UpSetR::upset(
      fromList(lst.protCond),
      nsets = length(unique(evidencekeys$condition)),
      nintersects = 50,
      order.by = "freq"
    )
    
    if(printPDF) garbage <- dev.off()
    
    if(verbose) message(" done ")
  }
  
  
  # Peptide ion oversampling
  if (plotPIO) {
    if(verbose) message("--- Plot Plot Ion Oversampling")
    
    if(printPDF) pdf(
      'QC_Plots_PepIonOversampling.pdf',
      width = nsamples * 3,
      height = 20,
      onefile = TRUE
    )
    fa <-
      ggplot(
        oversampling0[with(oversampling0, order(MSMS.counts)), ],
        aes(
          x = bioreplicate,
          y = FxOverSamp,
          fill = MSMS.counts,
          label = FxOverSamp
        )
      ) +
      geom_bar(stat = "identity", alpha = 0.7) +
      geom_text(
        aes(label = round(FxOverSamp, digits = 1)),
        vjust = 1 ,
        size = 10,
        position = position_stack()
      ) +
      xlab("Experiment") + ylab("Fraction (percentage)") +
      ggtitle("Peptide ion oversampling, all peptides reported by MaxQuant") +
      theme(legend.text = element_text(size = 20)) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 20
      )) +
      theme(axis.text.y = element_text(size = 20)) +
      theme(axis.title.x = element_text(size = 30)) +
      theme(axis.title.y = element_text(size = 30)) +
      theme(plot.title = element_text(size = 40)) +
      scale_fill_brewer(palette = "Set1")
    print(fa)
    
    fb <-
      ggplot(
        oversampling1[with(oversampling1, order(MSMS.counts)), ],
        aes(
          x = bioreplicate,
          y = FxOverSamp,
          fill = MSMS.counts,
          label = FxOverSamp
        )
      ) +
      geom_bar(stat = "identity", alpha = 0.7) +
      geom_text(
        aes(label = round(FxOverSamp, digits = 1)),
        vjust = 1 ,
        size = 10,
        position = position_stack()
      ) +
      xlab("Experiment") + ylab("Fraction (percentage)") +
      ggtitle("Peptide ion oversampling, only peptides detected (MS1 AUC calculated)") +
      theme(legend.text = element_text(size = 20)) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 20
      )) +
      theme(axis.text.y = element_text(size = 20)) +
      theme(axis.title.x = element_text(size = 30)) +
      theme(axis.title.y = element_text(size = 30)) +
      theme(plot.title = element_text(size = 40)) +
      scale_fill_brewer(palette = "Set1")
    print(fb)
    
    fc <-
      ggplot(
        oversampling2[with(oversampling2, order(MSMS.counts)), ],
        aes(
          x = bioreplicate,
          y = FxOverSamp,
          fill = MSMS.counts,
          label = FxOverSamp
        )
      ) +
      geom_bar(stat = "identity", alpha = 0.7) +
      geom_text(
        aes(label = round(FxOverSamp, digits = 1)),
        vjust = 1 ,
        size = 10,
        position = position_stack()
      ) +
      xlab("Experiment") + ylab("Fraction (percentage)") +
      ggtitle(
        "Peptide ion oversampling, only peptides detected (MS1 AUC calculated) and identified (confidence MS/MS)"
      ) +
      theme(legend.text = element_text(size = 20)) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 20
      )) +
      theme(axis.text.y = element_text(size = 20)) +
      theme(axis.title.x = element_text(size = 30)) +
      theme(axis.title.y = element_text(size = 30)) +
      theme(plot.title = element_text(size = 40)) +
      scale_fill_brewer(palette = "Set1")
    print(fc)
    
    if(printPDF) garbage <- dev.off()
    
    if(verbose) message(" done ")
  }
  
  
  # Charge state distribution
  if (plotCS) {
    if(verbose) message("--- Plot Charge State")
    
    if(printPDF) pdf(
      'QC_Plots_CHARGESTATE.pdf',
      width = nsamples * 3,
      height = 20,
      onefile = TRUE
    )
    ga <-
      ggplot(chargeState[with(chargeState, order(charge)), ],
             aes(x = bioreplicate, y = FxOverSamp, fill = charge)) +
      geom_bar(stat = "identity", alpha = 0.7) +
      ggrepel::geom_text_repel(
        aes(label = round(FxOverSamp, digits = 1)),
        vjust = -0.5,
        size = 5,
        position = position_stack()
      ) +
      xlab("Experiment") + ylab("Fraction (percentage)") +
      ggtitle("Precursor charge state distribution") +
      theme(legend.text = element_text(size = 20)) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 20
      )) +
      theme(axis.text.y = element_text(size = 20)) +
      theme(axis.title.x = element_text(size = 30)) +
      theme(axis.title.y = element_text(size = 30)) +
      theme(plot.title = element_text(size = 40)) +
      scale_fill_brewer(palette = "Spectral")
    print(ga)
    if(printPDF) garbage <- dev.off()
    
    if(verbose) message(" done ")
  }
  
  # Mass error
  if (plotME) {
    if(verbose) message("--- Plot Mass Error")
    if(printPDF) pdf(
      'QC_Plots_MASSERROR.pdf',
      width = nsamples * 3,
      height = 20,
      onefile = TRUE
    )
    fa <-
      ggplot(evidencekeys,
             aes(x = bioreplicate, y = uncalibrated.mass.error..ppm.)) +
      geom_boxplot(varwidth = TRUE, aes(fill = factor(condition)), alpha =
                     0.7) +
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
        size = 20
      ) +
      xlab("Experiment") + ylab("mass error") +
      ggtitle("Precursor mass error (in ppm) distribution, global median mass error on the top") +
      theme(legend.text = element_text(size = 20)) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 20
      )) +
      theme(axis.text.y = element_text(size = 20)) +
      theme(axis.title.x = element_text(size = 30)) +
      theme(axis.title.y = element_text(size = 30)) +
      theme(plot.title = element_text(size = 40)) +
      scale_fill_brewer(palette = "Spectral")
    print(fa)
    if(printPDF) garbage <- dev.off()
    if(verbose) message(" done ")
  }
  
  
  # mass-over-charge distribution
  if (plotMOCD) {
    if(verbose) message("--- Plot Mass-over-Charge distribution")
    
    if(printPDF) pdf(
      'QC_Plots_MZ.pdf',
      width = nsamples * 3,
      height = 20,
      onefile = TRUE
    )
    ga <- ggplot(evidencekeys, aes(x = bioreplicate, y = m.z)) +
      geom_boxplot(varwidth = TRUE, aes(fill = factor(condition)), alpha =
                     0.7) +
      geom_text(
        data = aggregate(m.z ~ bioreplicate, evidencekeys, median),
        aes(
          label = round(m.z, digits = 1),
          y = max(evidencekeys$m.z, na.rm = TRUE) + 30
        ),
        size = 20
      ) +
      xlab("Experiment") + ylab("m/z") +
      ggtitle("Precursor mass-over-charge distribution, global median m/z on the top") +
      theme(legend.text = element_text(size = 20)) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 20
      )) +
      theme(axis.text.y = element_text(size = 20)) +
      theme(axis.title.x = element_text(size = 30)) +
      theme(axis.title.y = element_text(size = 30)) +
      theme(plot.title = element_text(size = 40)) +
      scale_fill_brewer(palette = "Spectral")
    print(ga)
    if(printPDF) garbage <- dev.off()
    
    if(verbose) message(" done ")
  }
  
  
  #m Peptide Intensity CV
  if (plotPEPICV) {
    if(verbose) message("--- Plot Peptide Intensity CV")
    if(printPDF) pdf(
      'QC_Plots_PEPINT.pdf',
      width = nsamples * 3,
      height = 20,
      onefile = TRUE
    )
    peptCV <-
      data.table::data.table(subset(evidencekeys,!is.na(intensity)))
    peptCV <-
      peptCV[, list(intensity = sum(intensity, na.rm = TRUE)),
             by = list(condition, bioreplicate, modified.sequence)]
    peptCV <-
      peptCV[, list(
        pCV = 100 * (sd(intensity) / mean(intensity)),
        sumInt = sum(intensity),
        pDX = length(unique(bioreplicate))
      ),
      by = list(condition, modified.sequence)]
    peptCV <- peptCV[, bin.all := .artms_qcut(sumInt, 4)]
    peptCV <- peptCV[, bin.condition := .artms_qcut(sumInt, 4),
                     by = condition]
    
    ha <-
      ggplot(subset(peptCV,!is.na(pCV)), aes(condition, pCV)) +
      geom_boxplot(varwidth = TRUE, aes(fill = factor(condition)), alpha =
                     0.7) +
      geom_text(
        data = aggregate(pCV ~ condition, subset(peptCV,!is.na(pCV)), median),
        aes(
          label = round(pCV, digits = 1),
          y = max(peptCV$pCV, na.rm = TRUE) + 1
        ),
        size = 10
      ) +
      geom_text(
        data = aggregate(pCV ~ condition, subset(peptCV,!is.na(pCV)), length),
        aes(label = round(pCV, digits = 1), y = 0),
        size = 10
      ) +
      xlab("Condition") + ylab("Coefficient of variance (%)") +
      ggtitle(
        "Distribution of peptide feature intensity CV within each condition\n  
        Overall median CV for each condition is given on the top and number\n 
        of features used to calculate CVs is shown on the bottom"
      ) +
      theme(legend.text = element_text(size = 20)) +
      theme(axis.text.x = element_text(angle = 0, size = 20)) +
      theme(axis.text.y = element_text(size = 20)) +
      theme(axis.title.x = element_text(size = 30)) +
      theme(axis.title.y = element_text(size = 30)) +
      theme(plot.title = element_text(size = 40)) +
      scale_fill_brewer(palette = "Spectral")
    print(ha)
    
    hb <-
      ggplot(subset(peptCV,!is.na(pCV)), aes(interaction(condition, bin.condition), pCV)) +
      geom_boxplot(varwidth = TRUE, aes(fill = factor(condition)), alpha =
                     0.7) +
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
        size = 7,
        angle = 90
      ) +
      geom_text(
        data = aggregate(
          pCV ~ condition + bin.condition,
          subset(peptCV,!is.na(pCV)),
          length
        ),
        aes(label = round(pCV, digits = 1), y = 0),
        size = 7,
        angle = 90
      ) +
      xlab("Condition") + ylab("Coefficient of variance (%)") +
      ggtitle(
        "Distribution of peptide feature intensity CV within each condition\n  
        For each condition, peptides were ranked by summed intensity and the CV for each peptide was calculated,\n  
        therefore each condition shows 4 distribution (box) for low (1) to high (4) intensity peptides\n
        Overall median CV within each bin/condition is shown on the top and number of features used to calculate CV is given on the bottom"
      ) +
      theme(legend.text = element_text(size = 20)) +
      theme(axis.text.x = element_text(angle = 90, size = 20)) +
      theme(axis.text.y = element_text(size = 20)) +
      theme(axis.title.x = element_text(size = 30)) +
      theme(axis.title.y = element_text(size = 30)) +
      theme(plot.title = element_text(size = 20, hjust = 0.5)) +
      scale_fill_brewer(palette = "Spectral")
    print(hb)
    
    if(printPDF) garbage <- dev.off()
    if(verbose) message(" done ")
  }
  
  # peptide detection (using modified.sequence)
  if (plotPEPDETECT) {
    if(verbose) message("--- Plot Peptide Detection (using modified.sequence)")
    if(printPDF) pdf(
      'QC_Plots_PepDetect.pdf',
      width = nsamples * 3,
      height = 20,
      onefile = TRUE
    )
    peptDX <- peptCV[, .N, by = list(condition, pDX)]
    peptDXT <-
      peptDX[, list(N.total = sum(N)), by = list(condition)]
    peptDX <- merge(peptDX, peptDXT, by = "condition")
    peptDX$fxDx <- 100 * (peptDX$N / peptDX$N.total)
    
    ia <-
      ggplot(peptDX, aes(
        x = condition,
        y = fxDx,
        fill = factor(pDX)
      )) +
      geom_bar(stat = "identity", alpha = 0.7) +
      ggrepel::geom_text_repel(
        aes(label = round(fxDx, digits = 1)),
        vjust = 1,
        size = 10,
        position = position_stack()
      ) +
      xlab("Condition") + ylab("Counts") +
      ggtitle("Frequency of peptides detection across replicates within the same condition") +
      theme(legend.text = element_text(size = 20)) +
      theme(axis.text.x = element_text(angle = 0, size = 20)) +
      theme(axis.text.y = element_text(size = 20)) +
      theme(plot.title = element_text(size = 40)) +
      scale_fill_brewer(palette = "Set1")
    print(ia)
    
    if(printPDF) garbage <- dev.off()
    if(verbose) message(" done ")
  }
  
  
  #m Protein Intensity CV
  if (plotPROTICV) {
    if(verbose) message("--- Plot Protein Intensity CV")
    if(printPDF) pdf(
      'QC_Plots_ProtInt.pdf',
      width = nsamples * 3,
      height = 20,
      onefile = TRUE
    )
    protCV <-
      data.table::data.table(subset(evidencekeys,!is.na(intensity)))
    protCV <-
      protCV[, list(intensity = sum(intensity, na.rm = TRUE)), by = list(condition, bioreplicate, proteins)]
    protCV <-
      protCV[, list(
        pCV = 100 * (sd(intensity) / mean(intensity)),
        sumInt = sum(intensity),
        pDX = length(unique(bioreplicate))
      ), by = list(condition, proteins)]
    protCV <- protCV[, bin.all := .artms_qcut(sumInt, 4)]
    protCV <-
      protCV[, bin.condition := .artms_qcut(sumInt, 4), by = condition]
    
    ja <-
      ggplot(subset(protCV,!is.na(pCV)), aes(condition, pCV)) +
      geom_boxplot(varwidth = TRUE, aes(fill = factor(condition)), alpha =
                     0.7) +
      geom_text(
        data = aggregate(pCV ~ condition, subset(protCV,!is.na(pCV)), median),
        aes(
          label = round(pCV, digits = 1),
          y = max(protCV$pCV, na.rm = TRUE) + 1
        ),
        size = 15
      ) +
      xlab("Condition") + ylab("Coefficient of variance (%)") +
      ggtitle(
        "Distribution of Protein intensity CV within each condition\n
        Overall median CV for each condition is given on the top and\n
        number of proteins used to calculate CVs is shown on the bottom"
      ) +
      theme(legend.text = element_text(size = 20)) +
      theme(axis.text.x = element_text(angle = 0, size = 20)) +
      theme(axis.text.y = element_text(size = 20)) +
      theme(axis.title.x = element_text(size = 30)) +
      theme(axis.title.y = element_text(size = 30)) +
      theme(plot.title = element_text(size = 20, hjust = 0.5)) +
      scale_fill_brewer(palette = "Spectral")
    print(ja)
    
    jb <-
      ggplot(subset(protCV,!is.na(pCV)), aes(interaction(condition, bin.condition), pCV)) +
      geom_boxplot(varwidth = TRUE, aes(fill = factor(condition)), alpha =
                     0.7) +
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
        size = 7,
        angle = 90
      ) +
      geom_text(
        data = aggregate(
          pCV ~ condition + bin.condition,
          subset(protCV,!is.na(pCV)),
          length
        ),
        aes(label = round(pCV, digits = 1), y = 0),
        size = 7,
        angle = 90
      ) +
      xlab("Condition") + ylab("Coefficient of variance (%)") +
      ggtitle(
        "Distribution of Protein (summed) intensity CV within each condition\n 
        For each condition, Proteins were ranked by summed intensity and the CV for each protein was calculated\n
        therefore each condition shows 4 distribution (box) for low (1) to high (4) intensity proteins\n
        Overall median CV within each bin/condition is shown on the top and number of protein groups used to calculate CV is given on the bottom"
      ) +
      theme(legend.text = element_text(size = 20)) +
      theme(axis.text.x = element_text(angle = 90, size = 20)) +
      theme(axis.text.y = element_text(size = 20)) +
      theme(axis.title.x = element_text(size = 30)) +
      theme(axis.title.y = element_text(size = 30)) +
      theme(plot.title = element_text(size = 20, hjust = 0.5)) +
      scale_fill_brewer(palette = "Spectral")
    print(jb)
    
    if(printPDF) garbage <- dev.off()
    if(verbose) message(" done ")
  }
  
  # Protein detection
  if (plotPROTDETECT) {
    if(verbose) message("--- Plot Protein Detection")
    if(printPDF) pdf(
      'QC_Plots_ProtDetect.pdf',
      width = nsamples * 3,
      height = 20,
      onefile = TRUE
    )
    protDX <- protCV[, .N, by = list(condition, pDX)]
    protDXT <-
      protDX[, list(N.total = sum(N)), by = list(condition)]
    protDX <- merge(protDX, protDXT, by = "condition")
    protDX$fxDx <- 100 * (protDX$N / protDX$N.total)
    ka <-
      ggplot(protDX, aes(
        x = condition,
        y = fxDx,
        fill = factor(pDX)
      )) +
      geom_bar(stat = "identity", alpha = 0.7) +
      ggrepel::geom_text_repel(
        aes(label = round(fxDx, digits = 1)),
        vjust = 1,
        size = 10,
        position = position_stack()
      ) +
      xlab("Condition") + ylab("Counts") +
      ggtitle("Detection of Protein across replicates of the same condition") +
      theme(legend.text = element_text(size = 20)) +
      theme(axis.text.x = element_text(angle = 0, size = 20)) +
      theme(axis.text.y = element_text(size = 20)) +
      scale_fill_brewer(palette = "Set1")
    print(ka)
    
    kb <-
      ggplot(subset(evidencekeys,!is.na(intensity)), aes(bioreplicate, log2(intensity))) +
      geom_boxplot(varwidth = TRUE,
                   aes(fill = potential.contaminant),
                   alpha = 0.7) +
      #ggrepel::geom_text_repel(data = aggregate(intensity ~ bioreplicate + potential.contaminant, subset(evidencekeys, !is.na(intensity)), median), aes(label = round(log2(intensity), digits=1), y = log2(max(evidencekeys$intensity, na.rm= TRUE))+0.5 ), size = 15) +
      xlab("Experiment") + ylab("Log2 Intensity") +
      ggtitle("Peptide feature intensity distribution") +
      theme(legend.text = element_text(size = 20)) +
      theme(axis.text.x = element_text(angle = 90, size = 20)) +
      theme(axis.text.y = element_text(size = 20)) +
      theme(axis.title.x = element_text(size = 30)) +
      theme(axis.title.y = element_text(size = 30)) +
      theme(plot.title = element_text(size = 40)) +
      scale_fill_brewer(palette = "Set1")
    print(kb)
    
    if ("fraction" %in% colnames(evidencekeys)) {
      kc <-
        ggplot(subset(evidencekeys,!is.na(intensity)), aes(bioreplicate, log2(intensity))) +
        geom_boxplot(varwidth = TRUE,
                     aes(fill = potential.contaminant),
                     alpha = 0.7) +
        #ggrepel::geom_text_repel(data = aggregate(intensity ~ bioreplicate + potential.contaminant, subset(evidencekeys, !is.na(intensity)), median), aes(label = round(log2(intensity), digits=1), y = log2(max(evidencekeys$intensity, na.rm= TRUE))+0.5 ), size = 15) +
        xlab("Experiment") + ylab("Log2 Intensity") +
        facet_wrap( ~ fraction, ncol = 5) +
        ggtitle("Peptide feature intensity distribution by Fraction") +
        theme(legend.text = element_text(size = 20)) +
        theme(axis.text.x = element_text(angle = 90, size = 20)) +
        theme(axis.text.y = element_text(size = 20)) +
        theme(axis.title.x = element_text(size = 30)) +
        theme(axis.title.y = element_text(size = 30)) +
        theme(plot.title = element_text(size = 40)) +
        scale_fill_brewer(palette = "Set1")
      print(kc)
      
    }
    
    kd <- ggplot(evidencekeys, aes(x = log2(intensity))) +
      geom_density(alpha = .5, aes(fill = type)) +
      ggtitle("Peptide feature intensity distribution by ID Type") +
      facet_wrap( ~ bioreplicate, ncol = 5) +
      theme(legend.text = element_text(size = 20)) +
      theme(axis.text.x = element_text(size = 20)) +
      theme(axis.text.y = element_text(size = 20)) +
      theme(axis.title.x = element_text(size = 30)) +
      theme(axis.title.y = element_text(size = 30)) +
      theme(plot.title = element_text(size = 40)) +
      scale_fill_brewer(palette = "Set1")
    print(kd)
    if(printPDF) garbage <- dev.off()
    if(verbose) message(" done ")
  }
  
  
  if (plotIDoverlap) {
    if(verbose) message("--- Plot ID overlap")
    if(printPDF) pdf(
      'QC-ID-Overlap.pdf',
      width = nsamples * 4,
      height = nsamples * 4,
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
    
    gplots::heatmap.2(
      ovlSeqM * 100,
      col = hmcol,
      Rowv = FALSE,
      Colv = "Rowv",
      cexRow = 5,
      cexCol = 5,
      margins = c(40, 40),
      cellnote = round(ovlSeqM * 100, 1),
      notecex = 8,
      main = "Pairwise peptide identification overlap (only peptides with at least 1 peptide detected and identified)",
      notecol = "black",
      trace = "none",
      key = TRUE,
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
      dendrogram = "none"
    )
    
    if ("fraction" %in% colnames(evidencekeys)) {
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
        main = "Pairwise peptide identification overlap by Fraction (only peptides with at least 1 peptide detected and identified)",
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
      cexRow = 5,
      cexCol = 5,
      margins = c(40, 40),
      cellnote = round(ovlProtM * 100, 1),
      notecex = 8,
      main = "Pairwise protein identification overlap (only proteins with at least 1 peptide detected and identified)",
      notecol = "black",
      trace = "none",
      key = TRUE,
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
  
  if (plotIC) {
    if(verbose) message("--- Plot Inter-Correlation")
    if(printPDF) pdf(
      'QC-IntCorrelation.pdf',
      width = nsamples * 4,
      height = nsamples * 4,
      onefile = TRUE
    )
    peptintmtx <-
      subset(
        evidencekeys,
        !is.na(intensity),
        select = c("bioreplicate", "sequence", "intensity")
      )
    peptintmtx <-
      data.table::dcast(peptintmtx, sequence ~ bioreplicate, sum, value.var = "intensity")
    peptintmtx <- as.matrix(peptintmtx[, -1])
    peptintmtx <- log2(peptintmtx)
    peptintmtx[!is.finite(peptintmtx)] <- NA
    pairs(
      peptintmtx,
      upper.panel = .artms_panelCor,
      diag.panel = .artms_panelHist,
      lower.panel = panel.smooth,
      pch = "."
    )
    
    mpept <- peptintmtx
    mpept[is.na(mpept)] <- 0
    pept.pca <- prcomp(t(mpept))
    pept.pcax <- as.data.frame(pept.pca$x)
    pept.pcax <-
      merge(unique(evidencekeys[, c("bioreplicate", "condition")]),
            pept.pcax,
            by.x = "bioreplicate",
            by.y = "row.names")
    
    ma <-
      ggplot(pept.pcax, aes(PC1, PC2)) + #, color=condition, shape=condition
      geom_point(alpha = .8, size = 30) +
      geom_label_repel(
        aes(label = bioreplicate),
        #label.size = 10,
        box.padding   = 0.1,
        point.padding = 0.1,
        segment.color = 'grey50'
      ) +
      ggtitle("Peptide Intensity Principal Component Analysis")
    # ggalt is giving problem: it cannot be installed on linux
    # ggalt::geom_encircle(aes(group=condition, fill=condition),alpha=0.4)
    print(ma)
    
    protintmtx <-
      subset(
        evidencekeys,
        !is.na(intensity),
        select = c("bioreplicate", "proteins", "intensity")
      )
    protintmtx <-
      data.table::dcast(protintmtx, proteins ~ bioreplicate, sum, value.var = "intensity")
    protintmtx <- as.matrix(protintmtx[, -1])
    protintmtx <- log2(protintmtx)
    protintmtx[!is.finite(protintmtx)] <- NA
    pairs(
      protintmtx,
      upper.panel = .artms_panelCor,
      diag.panel = .artms_panelHist,
      lower.panel = panel.smooth,
      pch = "."
    )
    
    mprot <- protintmtx
    mprot[is.na(mprot)] <- 0
    prot.pca <- prcomp(t(mprot))
    prot.pcax <- as.data.frame(prot.pca$x)
    prot.pcax <-
      merge(unique(evidencekeys[, c("bioreplicate", "condition")]),
            prot.pcax,
            by.x = "bioreplicate",
            by.y = "row.names")
    
    mb <-
      ggplot(prot.pcax, aes(PC1, PC2)) + #, color=condition, shape=condition
      geom_point(alpha = .8, size = 40) +
      geom_label_repel(
        aes(label = bioreplicate),
        #label.size = 10,
        box.padding   = 0.1,
        point.padding = 0.1,
        segment.color = 'grey50'
      ) +
      ggtitle("Protein Intensity Principal Component Analysis")
    # ggalt is giving problem: it cannot be installed on linux
    # ggalt::geom_encircle(aes(group=condition, fill=condition),alpha=0.4)
    print(mb)
    
    if(printPDF) garbage <- dev.off()
    if(verbose) message(" done ")
  }
  
  
  if (plotSP) {
    if(verbose) message("--- Plot Sample Preparation")
    if(printPDF) pdf(
      'QC-SamplePrep.pdf',
      width = nsamples * 4,
      height = 20,
      onefile = TRUE
    )
    evidence.misscleavages <-
      data.table(subset(evidencekeys, ms.ms.count > 0))
    
    misscleavages.tot <-
      evidence.misscleavages[, .N, by = list(bioreplicate)]
    
    misscleavages.dt <- evidence.misscleavages[, .N,
                                               by = list(bioreplicate,
                                                         missed.cleavages)]
    misscleavages.dt <- merge(misscleavages.dt,
                              misscleavages.tot,
                              by = "bioreplicate",
                              all = TRUE)
    
    misscleavages.dt$mc <-
      as.numeric(format(
        100 * (misscleavages.dt$N.x / misscleavages.dt$N.y),
        digits = 5
      ))
    
    na <-
      ggplot(misscleavages.dt,
             aes(
               x = bioreplicate,
               y = mc,
               fill = as.factor(missed.cleavages),
               label = mc
             )) +
      geom_bar(stat = "identity", alpha = 0.7) +
      geom_text(
        aes(label = round(mc, digits = 1)),
        vjust = 0,
        size = 10,
        position = position_stack()
      ) +
      xlab("Experiment") + ylab("Fraction (percentage)") +
      ggtitle("Missing cleavage stats") +
      theme(legend.text = element_text(size = 20)) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 20
      )) +
      theme(axis.text.y = element_text(size = 20)) +
      theme(axis.title.x = element_text(size = 30)) +
      theme(axis.title.y = element_text(size = 30)) +
      theme(plot.title = element_text(size = 40)) +
      scale_fill_brewer(palette = "Set1")
    print(na)
    
    
    oxMet.df <-
      plyr::ddply(
        evidencekeys,
        c("bioreplicate", "condition"),
        plyr::summarise,
        pct.OxM = (length(oxidation..m.[oxidation..m. > 0]) / length(sequence)) *
          100
      )
    
    nb <- ggplot(oxMet.df,
                 aes(
                   x = bioreplicate,
                   y = pct.OxM,
                   fill = condition,
                   label = pct.OxM
                 )) +
      geom_bar(stat = "identity", alpha = 0.7) +
      geom_text(
        aes(label = round(pct.OxM,
                          digits = 1)),
        vjust = 0,
        size = 15,
        position = position_stack()
      ) +
      xlab("Experiment") +
      ylab("Fraction (percentage)") +
      ggtitle("Percentage of peptides where at least 1 Methionine is oxidized") +
      theme(legend.text = element_text(size = 20)) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 20
      )) +
      theme(axis.text.y = element_text(size = 20)) +
      theme(axis.title.x = element_text(size = 30)) +
      theme(axis.title.y = element_text(size = 30)) +
      theme(plot.title = element_text(size = 40)) +
      scale_fill_brewer(palette = "Spectral")
    print(nb)
    if(printPDF) garbage <- dev.off()
    if(verbose) message(" done ")
  }
  
  
  
} # END OF artmsQualityControlEvidenceExtended  


# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
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
  m <-
    matrix(
      unlist(pepmatch),
      nrow = length(pept),
      ncol = length(pept),
      byrow = TRUE
    )
  colnames(m) <- rownames(m) <- names(pept)
  ans <- list(Peptides = pept, M = m)
  return(ans)
}

# ------------------------------------------------------------------------------
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
.artms_panelCor <-
  function(x,
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

# ------------------------------------------------------------------------------
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