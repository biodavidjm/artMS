# ------------------------------------------------------------------------------
# @title Run MSstats
# @description Run MSstats giving a processed evidence file (MSstats format)
# and a contrast file. It also generates a series of summary plots before,
# and after normalization.
# @param dmss (data.frame) Formatted and filtered evidence file 
# (MSstats format)
# @param contrasts (data.frame) The contrast data.frame in MSstats format
# @param config (yaml object) the configation (imported yaml) object
# @param printPDF (logical) if `TRUE` (default) prints pdf files
# @param printTables (logical) `TRUE` (default) prints results tables
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
                              printPDF = TRUE,
                              printTables = TRUE,
                              verbose = TRUE) {
  
  # Check configuration parameters for msstats-----
  config <- .artms_provide_msstats_config_miss_parameters(config = config, 
                                                          verbose = verbose)

  # plot the data BEFORE normalization-----
  if(printPDF){
    if (grepl('before', config$msstats$profilePlots)) {
      if(verbose) message("-- QC PLOT: before")
      mssquant <- dataProcess(
        raw = dmss,
        logTrans = config$msstats$logTrans,
        normalization = FALSE,
        nameStandards = NULL,
        featureSubset = config$msstats$feature_subset,
        remove_uninformative_feature_outlier = config$msstats$remove_uninformative_feature_outlier,
        min_feature_count = config$msstats$min_feature_count,
        n_top_feature = config$msstats$n_top_feature,
        summaryMethod = config$msstats$summaryMethod,
        equalFeatureVar = config$msstats$equalFeatureVar,
        censoredInt = config$msstats$censoredInt,
        MBimpute = config$msstats$MBimpute,
        remove50missing = config$msstats$remove50missing,
        fix_missing = unlist(config$msstats$fix_missing),
        maxQuantileforCensored = config$msstats$maxQuantileforCensored,
        use_log_file = config$msstats$use_log_file,
        append = config$msstats$append,
        verbose = verbose,
        log_file_path = unlist(config$msstats$log_file_path)
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
  } #if print pdf
  
  # Normalization------
  
  if (!is.null(config$msstats$normalization_reference) &
      config$msstats$normalization_method == 'globalStandards') {
    # if globalStandars is selected, must have a reference protein(s)
    normalization_refs = unlist(lapply(
      strsplit(config$msstats$normalization_reference, split = ','),
      FUN = .artms_trim
    ))
    mssquant <- dataProcess(
      raw = dmss,
      logTrans = config$msstats$logTrans,
      normalization = config$msstats$normalization_method,
      nameStandards = normalization_refs,
      featureSubset = config$msstats$feature_subset,
      remove_uninformative_feature_outlier = config$msstats$remove_uninformative_feature_outlier,
      min_feature_count = config$msstats$min_feature_count,
      n_top_feature = config$msstats$n_top_feature,
      summaryMethod = config$msstats$summaryMethod,
      equalFeatureVar = config$msstats$equalFeatureVar,
      censoredInt = config$msstats$censoredInt,
      MBimpute = config$msstats$MBimpute,
      remove50missing = config$msstats$remove50missing,
      fix_missing = unlist(config$msstats$fix_missing),
      maxQuantileforCensored = config$msstats$maxQuantileforCensored,
      use_log_file = config$msstats$use_log_file,
      append = config$msstats$append,
      verbose = verbose,
      log_file_path = unlist(config$msstats$log_file_path)
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
      logTrans = config$msstats$logTrans,
      normalization = config$msstats$normalization_method,
      nameStandards = NULL,
      featureSubset = config$msstats$feature_subset,
      remove_uninformative_feature_outlier = config$msstats$remove_uninformative_feature_outlier,
      min_feature_count = config$msstats$min_feature_count,
      n_top_feature = config$msstats$n_top_feature,
      summaryMethod = config$msstats$summaryMethod,
      equalFeatureVar = config$msstats$equalFeatureVar,
      censoredInt = config$msstats$censoredInt,
      MBimpute = config$msstats$MBimpute,
      remove50missing = config$msstats$remove50missing,
      fix_missing = unlist(config$msstats$fix_missing),
      maxQuantileforCensored = config$msstats$maxQuantileforCensored,
      use_log_file = config$msstats$use_log_file,
      append = config$msstats$append,
      verbose = verbose,
      log_file_path = unlist(config$msstats$log_file_path)
    )
  }
  
  # plot the data AFTER normalization-----
  if(printPDF){
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
  }
  
  if (!all(levels(mssquant$ProteinLevelData$GROUP) == colnames(contrasts))) {
    stop(
      sprintf(
        '\tERROR IN CONTRAST COMPARISON: GROUP LEVELS DIFFERENT FROM 
        CONTRASTS FILE \tGROUP LEVELS\t%s\n\tCONTRASTS FILE\t%s ',
        paste(levels(mssquant$GROUP_ORIGINAL), collapse = ','),
        paste(colnames(contrasts), collapse = ',')
      )
    )
  }
  
  # protein sample/group quantification----
  if(printTables){
    write.table(
      quantification(mssquant, type = "Sample"),
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
      mssquant$FeatureLevelData,
      file = gsub('.txt', '-mss-FeatureLevelData.txt', config$files$output),
      eol = "\n",
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE
    )
    
    write.table(
      mssquant$ProteinLevelData,
      file = gsub('.txt', '-mss-ProteinLevelData.txt', config$files$output),
      eol = "\n",
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE
    )
  }
  
  if(verbose) message("-- GROUP COMPARISON")
  results <- groupComparison(data = mssquant, 
                             contrast.matrix = contrasts, 
                             use_log_file = FALSE)
  
  if(printTables){
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
  }

  if(verbose) message("-- MSstats FINISHED! ")
  
  #(1) Minimal number of biological replicates per condition-----
  if(verbose) message(">> CALCULATING SAMPLE SIZE FOR FUTURE EXPERIMENTS ")
  results.ss1 <- designSampleSize(data = results$FittedModel,
                                  numSample = TRUE,
                                  desiredFC = c(0.58, 2),
                                  FDR = 0.05,
                                  power = 0.95, 
                                  use_log_file = FALSE)
    
    
  results.ss2 <- designSampleSize(data = results$FittedModel,
                                  numSample = TRUE,
                                  desiredFC = c(0.58, 2),
                                  FDR = 0.05,
                                  power = 0.9, 
                                  use_log_file = FALSE)
  
  results.ss3 <- designSampleSize(data = results$FittedModel,
                                  numSample = TRUE,
                                  desiredFC = c(0.58, 2),
                                  FDR = 0.05,
                                  power = 0.8, 
                                  use_log_file = FALSE)
    
    
  results.ss <- rbind(results.ss1, results.ss2, results.ss3)

  #(2) Power calculation-----
  if(verbose) message(">> CALCULATING POWER OF EXPERIMENT")
  results.power1 <- designSampleSize(data = results$FittedModel,
                                     numSample = 2,
                                     desiredFC = c(0.58, 2),
                                     FDR = 0.05,
                                     power = TRUE, 
                                     use_log_file = FALSE)
  results.power2 <- designSampleSize(data = results$FittedModel,
                                     numSample = 3,
                                     desiredFC = c(0.58, 2),
                                     FDR = 0.05,
                                     power = TRUE, 
                                     use_log_file = FALSE)
  results.power3 <- designSampleSize(data = results$FittedModel,
                                     numSample = 4,
                                     desiredFC = c(0.58, 2),
                                     FDR = 0.05,
                                     power = TRUE, 
                                     use_log_file = FALSE)
  results.power4 <- designSampleSize(data = results$FittedModel,
                                     numSample = 5,
                                     desiredFC = c(0.58, 2),
                                     FDR = 0.05,
                                     power = TRUE, 
                                     use_log_file = FALSE)

  results.power <- rbind(results.power1, results.power2, results.power3, results.power4)
  
  if(printTables){
    write.table(
      results.ss,
      file = gsub(".txt", "_sampleSize.txt", config$files$output),
      eol = "\n",
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE
    )
    write.table(
      results.power,
      file = gsub(".txt", "_experimentPower.txt", config$files$output),
      eol = "\n",
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE
    )
  }

  results$power <- results.power
  results$sample_size <- results.ss
  
  return(results)
}