# ------------------------------------------------------------------------------
# @title Write extras
# @description Extras after MSstats, as annotations, volcano plots, heatmaps
# @param results MSstats results
# @param config The configuration object (yaml)
# @return Extras as selected in the yaml file, including:
# - volcano plot (pdf)
# - Adding annotations (gene symbol based on uniprot)
# @keywords extras, annotations, volcano
.artms_writeExtras <- function(results, 
                               config) {
  
  if (length(results) == 0 | !exists('results')) {
    stop("ERROR!! NO RESULTS FOUND TO ANNOTATE!")
  }
  
  # Annotation
  if (config$output_extras$annotate$enabled) {
    results_ann <- artms_annotationUniprot(x = results, 
                                columnid = 'Protein', 
                                species = config$output_extras$annotate$species)
    output_annotated_file <- gsub(".txt", "-annotated.txt", config$files$output)
    write.table(results_ann, output_annotated_file, quote = FALSE, 
                row.names = FALSE, 
                col.names = TRUE, 
                sep = "\t")
  }else{
    results_ann <- results
  }
  
  lfc_lower <-
    as.numeric(unlist(strsplit(config$output_extras$plots$LFC, split = " "))[1])
  lfc_upper <-
    as.numeric(unlist(strsplit(config$output_extras$plots$LFC, split = " "))[2])
  
  ## This option was originally available in the config file but
  ## it was removed for now
  config$output_extras$plots$comparisons <- "all"
  ## ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  
  selected_labels <- config$output_extras$plots$comparisons
  if (is.null(selected_labels) ||
      selected_labels == 'all')
    selected_labels = '*'
  
  # remove the Inf, -Inf log2FC hits.
  results_ann <- results_ann[!is.infinite(results_ann$log2FC), ]
  
  ## select data points  by LFC & FDR criterium in single condition and
  ## adding corresponding data points from the other conditions
  sign_hits <-
    .artms_significantHits(
      results_ann,
      labels = selected_labels,
      LFC = c(lfc_lower, lfc_upper),
      FDR = config$output_extras$plots$FDR
    )
  if (dim(sign_hits)[1] == 0)
    stop("NO SIGNIFICANT HITS DETECTED IN THIS EXPERIMENT. ABORTING PLOTS.\n")
  sign_labels <- unique(sign_hits$Label)
  cat(
    sprintf(
      "\tSELECTED HITS FOR PLOTS WITH LFC BETWEEN %s AND %s AT %s FDR:\t%s\n",
      lfc_lower,
      lfc_upper,
      config$output_extras$plots$FDR,
      nrow(sign_hits) / length(sign_labels)
    )
  )
  
  ## REPRESENTING RESULTS AS HEATMAP, if enabled
  if (config$output_extras$plots$heatmap) {
    # Heatmap only for > 1 comparison
    if (dim(sign_hits)[1] > 1) {
      cat(">> PLOTTING HEATMAP FOR SIGNIFICANT CHANGES\n")
      heat_labels <-
        .artms_prettyPrintHeatmapLabels(
          uniprot_acs = sign_hits$Protein,
          uniprot_ids = sign_hits$name,
          gene_names = sign_hits$Gene.names
        )
      heat_data_w <-
        .artms_plotHeat(
          mss_F = sign_hits,
          out_file =  gsub('.txt', '-sign.pdf', config$files$output),
          names = heat_labels,
          cluster_cols = config$output_extras$plots$heatmap_cluster_cols,
          display = config$output_extras$plots$heatmap_display
        )
    }
  }
  
  if (config$output_extras$plots$volcano) {
    cat(">> PLOTTING VOLCANO PLOT\n")
    file_name <- gsub('.txt', '-volcano.pdf', config$files$output)
    artms_volcanoPlot(
      mss_results = results_ann[grep(selected_labels, results_ann$Label), ], 
      lfc_upper = lfc_upper, lfc_lower = lfc_lower, 
      FDR = config$output_extras$plots$FDR,
      output_name = file_name
    )
  }
}



