# ------------------------------------------------------------------------------
# @title Run MSstats
# @description Run MSstats giving a processed evidence file (MSstats format)
# and a contrast file. It also generates a series of summary plots before,
# and after normalization.
# @param dmss (data.frame) Formatted and filtered evidence file 
# (MSstats format)
# @param contrasts (data.frame) The contrast data.frame in MSstats format
# @param config (yaml object) the configation (imported yaml) object
# @param verbose (logical) `TRUE` (default) shows function messages
# @return It generates several output files
# - If selected `before` and/or `after`, the `ProfilePlot` and `QCPlot` plots
# by the MSstats `dataProcessPlots` function are generated
# (in `.pdf` format)
# - Text file output of `quantification()` (`-mss-sampleQuant.txt`)
# - Text file output of `quantification(type="Group")` (`-mss-groupQuant.txt`)
# - MSstats `ProcessedData` normalized results (`-mss-normalized.txt`)
# - MSstats `ComparisonResult` results
# - MSstats `ModelQC` results
# - MSstats `designSampleSize` sample size
# - MSstats `designSampleSize` power experiment
# @keywords internal, run, MSstats, contrast, intensity, plots, QC
.artms_runMSstats <- function(dmss, 
                              contrasts, 
                              config,
                              verbose = TRUE) {
  
  # plot the data BEFORE normalization
  if (grepl('before', config$msstats$profilePlots)) {
    if(verbose) message("-- QC PLOT: before")
    mssquant <- dataProcess(
      raw = dmss,
      normalization = FALSE,
      #Interference has been depracated, but it was never used anyway
      # betweenRunInterferenceScore = config$msstats$interference,
      fillIncompleteRows = TRUE,
      summaryMethod = config$msstats$summaryMethod,
      censoredInt = config$msstats$censoredInt,
      cutoffCensored = config$msstats$cutoffCensored,
      MBimpute = config$msstats$MBimpute,
      featureSubset = config$msstats$feature_subset
    )
    dataProcessPlots(
      data = mssquant,
      type = "ProfilePlot",
      featureName = "Peptide",
      address = gsub('.txt', '-before', config$files$output)
    )
    dataProcessPlots(
      data = mssquant,
      type = "QCPlot",
      address = gsub('.txt', '-before', config$files$output)
    )
  }
  
  # Normalization
  
  if (!is.null(config$msstats$normalization_reference) &
      config$msstats$normalization_method == 'globalStandards') {
    # if globalStandars is selected, must have a reference protein(s)
    normalization_refs = unlist(lapply(
      strsplit(config$msstats$normalization_reference, split = ','),
      FUN = .artms_trim
    ))
    mssquant <- dataProcess(
      raw = dmss,
      normalization = config$msstats$normalization_method,
      nameStandards = normalization_refs,
      # Interference has been depracated, but it was never used anyway
      # betweenRunInterferenceScore = config$msstats$interference,
      fillIncompleteRows = TRUE,
      summaryMethod = config$msstats$summaryMethod,
      censoredInt = config$msstats$censoredInt,
      cutoffCensored = config$msstats$cutoffCensored,
      MBimpute = config$msstats$MBimpute,
      featureSubset = config$msstats$feature_subset
    )
  } else{
    if(verbose) message(
      sprintf(
        '-- Normalization method: %s',
        config$msstats$normalization_method
      )
    )
  
    mssquant = dataProcess(
      raw = dmss,
      normalization = config$msstats$normalization_method,
      #Interference has been depracated, but it was never used anyway
      # betweenRunInterferenceScore = config$msstats$interference,
      fillIncompleteRows = TRUE,
      summaryMethod = config$msstats$summaryMethod,
      censoredInt = config$msstats$censoredInt,
      cutoffCensored = config$msstats$cutoffCensored,
      MBimpute = config$msstats$MBimpute,
      featureSubset = config$msstats$feature_subset
    )
  }
  
  # plot the data AFTER normalization
  if (grepl('after', config$msstats$profilePlots)) {
    if(verbose) message("-- QC PLOT: after")
    dataProcessPlots(
      data = mssquant,
      type = "ProfilePlot",
      featureName = "Peptide",
      address = gsub('.txt', '-after', config$files$output)
    )
    dataProcessPlots(
      data = mssquant,
      type = "QCPlot",
      address = gsub('.txt', '-after', config$files$output)
    )
  }
  
  if (!all(levels(mssquant$GROUP_ORIGINAL) == colnames(contrasts))) {
    stop(
      sprintf(
        '\tERROR IN CONTRAST COMPARISON: GROUP LEVELS DIFFERENT FROM 
        CONTRASTS FILE \tGROUP LEVELS\t%s\n\tCONTRASTS FILE\t%s ',
        paste(levels(mssquant$GROUP_ORIGINAL), collapse = ','),
        paste(colnames(contrasts), collapse = ',')
      )
    )
  }
  
  # protein sample/group quantification
  write.table(
    quantification(mssquant),
    file = gsub('.txt', '-mss-sampleQuant.txt', config$files$output),
    eol = "\n",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  write.table(
    quantification(mssquant, type = "Group"),
    file = gsub('.txt', '-mss-groupQuant.txt', config$files$output),
    eol = "\n",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  
  if(verbose) message(sprintf('-- FITTING CONTRASTS:\t%s\n',paste(rownames(contrasts), collapse = ',')))
  
  write.table(
    mssquant$ProcessedData,
    file = gsub('.txt', '-mss-normalized.txt', config$files$output),
    eol = "\n",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  
  if(verbose) message("-- GROUP COMPARISON")
  results <- groupComparison(data = mssquant, contrast.matrix = contrasts)
  
  write.table(
    results$ComparisonResult,
    file = config$files$output,
    eol = "\n",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  write.table(
    results$ModelQC,
    file = gsub(".txt", "_ModelQC.txt", config$files$output),
    eol = "\n",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  write.table(
    mssquant$RunlevelData,
    file = gsub(".txt", "_RunLevelData.txt", config$files$output),
    eol = "\n",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )

  if(verbose) message("-- MSstats FINISHED! ")
  
  #(1) Minimal number of biological replicates per condition
  if(verbose) message(">> CALCULATING SAMPLE SIZE FOR FUTURE EXPERIMENTS ")
  results.ss1 <- designSampleSize(data = results$fittedmodel,
                                  numSample = TRUE,
                                  desiredFC = c(1.25, 2),
                                  FDR = 0.05,
                                  power = 0.95)
    
    
  results.ss2 <- designSampleSize(data = results$fittedmodel,
                                  numSample = TRUE,
                                  desiredFC = c(1.25, 2),
                                  FDR = 0.05,
                                  power = 0.9)
  
  results.ss3 <- designSampleSize(data = results$fittedmodel,
                                  numSample = TRUE,
                                  desiredFC = c(1.25, 2),
                                  FDR = 0.05,
                                  power = 0.8)
    
    
  results.ss <- rbind(results.ss1, results.ss2, results.ss3)
  write.table(
    results.ss,
    file = gsub(".txt", "_sampleSize.txt", config$files$output),
    eol = "\n",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  
  #(2) Power calculation
  if(verbose) message(">> CALCULATING POWER OF EXPERIMENT")
  results.power1 <-
    designSampleSize(
      data = results$fittedmodel,
      numSample = 3,
      desiredFC = c(1.25, 2),
      FDR = 0.05,
      power = TRUE
    )
  results.power2 <-
    designSampleSize(
      data = results$fittedmodel,
      numSample = 2,
      desiredFC = c(1.25, 2),
      FDR = 0.05,
      power = TRUE
    )
  results.power <- rbind(results.power1, results.power2)
  write.table(
    results.power,
    file = gsub(".txt", "_experimentPower.txt", config$files$output),
    eol = "\n",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  return(results)
}