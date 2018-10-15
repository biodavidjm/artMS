# ==============================================================================
## Small MSstats-related Functions

# ------------------------------------------------------------------------------
# @title Long to Wide format using the `Sequence` column of the evidence file
#
# @description Facilitates applying the dcast function, i.e., takes long-format
# data and casts it into wide-format data.
# @param d_long (data.frame) in long format
# @return (data.frame) Evidence file reshaped by rawfile and IsotopeLabelType
# @keywords internal, data.frame, dcast
.artms_castMaxQToWide <- function(d_long) {
  data_w <-
    data.table::dcast(
      Proteins + Sequence + Charge ~ RawFile + IsotopeLabelType,
      data = d_long,
      value.var = 'Intensity',
      fun.aggregate = sum,
      fill = NA
    )
  return(data_w)
}

# ------------------------------------------------------------------------------
# @title Long to Wide format selecting the `Modified.sequence` column of the
# evidence file
#
# @description Facilitates applying the dcast function, i.e., takes long-format
# data and casts it into wide-format data.
# @param d_long (data.frame) in long format
# @return (data.frame) Evidence file reshaped by rawfile and IsotopeLabelType
# @keywords internal, data.frame, dcast, ptm
.artms_castMaxQToWidePTM <- function(d_long) {
  data_w <-
    data.table::dcast(
      Proteins + Modified.sequence + Charge ~ RawFile + IsotopeLabelType,
      data = d_long,
      value.var = 'Intensity',
      fun.aggregate = sum,
      fill = NA
    )
  setnames(data_w, 2, 'Sequence')
  return(data_w)
}

# ------------------------------------------------------------------------------
# @title Check the `Raw file` column name on the evidence or keys data.frame
#
# @description Depending on how the data is loaded, the `Raw file` column
# might have different format. This function check to ensure consistency in
# both the evidence and keys data.frames
# @param (data.frame) keys or evidence files
# @return (data.frame) with the `RawFile` column name
# @keywords internal rawfile, columname
.artms_checkRawFileColumnName <- function(df) {
  if (!('RawFile' %in% colnames(df))) {
    if ("Raw.file" %in% colnames(df)) {
      df <- artms_changeColumnName(df, 'Raw.file', 'RawFile')
    } else if ("Raw file" %in% colnames(df)) {
      df <- artms_changeColumnName(df, 'Raw file', 'RawFile')
    } else{
      stop("\tERROR: CANNOT FIND THE Raw.file COLUMN\nPlease, revise it")
    }
  }
  return(df)
}

# ------------------------------------------------------------------------------
#' @title Change a specific column name in a given data.frame
#'
#' @description Making easier to change a column name in any data.frame
#' @param dataset (data.frame) with the column name you want to change
#' @param oldname (char) the old column name
#' @param newname (char) the new name for that column
#' @return (data.frame) with the new specified column name
#' @keywords rename, data.frame, columns
#' @examples
#' artms_data_ph_evidence <- artms_changeColumnName(
#'                                dataset = artms_data_ph_evidence,
#'                                oldname = "Phospho..STY.",
#'                                newname = "PH_STY")
#' @export
artms_changeColumnName <- function(dataset, oldname, newname) {
  if (!(oldname %in% colnames(dataset))) {
    stop("\nThe Column name provided <",
         oldname,
         "> was not found in the object provided\n")
  }
  colnames(dataset)[grep(paste0('^', oldname, '$'), colnames(dataset))] <-
    newname
  return(dataset)
}

# ------------------------------------------------------------------------------
# @title Filtering data
#
# @description Apply the filtering options, i.e., remove protein groups and/or
# contaminants, and/or, select posttranslational modification (if any)
# @param data (data.frame) Evidence file
# @param config (yaml.object) Configuration object (opened yaml file)
# @return (data.frame) filtered according to the options selected
# @keywords internal, filtering, remove, proteingroups, ptms
.artms_filterData <- function(data, config) {
  cat("\n>> FILTERING\n")
  if (config$data$filters$protein_groups == 'remove') {
    cat("\tPROTEIN GROUPS\tREMOVE\n")
    data_f <- .artms_removeMaxQProteinGroups(data)
  } else if (config$data$filters$protein_groups == 'keep') {
    cat("\tPROTEIN GROUPS\tIGNORE\n")
    data_f <- data
  } else{
    stop(
      "\n\nFILTERING OPTION FOR protein_groups 
      NOT UNDERSTOOD (OPTIONS AVAILABLE: keep OR remove\n\n"
    )
  }
  
  if (config$data$filters$contaminants) {
    cat("\tCONTAMINANTS\tREMOVE\n")
    data_f <- artms_filterEvidenceContaminants(data_f)
  }
  
  # DEAL WITH OLD CONFIGURATION FILES WHEN config$data$filters$modification 
  # COULD BE EMPTY
  if (is.null(config$data$filters$modification)) {
    cat("\tNO config$data$filters$modification provided. 
        Using 'AB' as default\n")
  } else if (config$data$filters$modification == 'AB' |
             config$data$filters$modification == 'APMS') {
    cat(sprintf("\tPROCESSING\t%s\n", config$data$filters$modification))
  } else if (config$data$filters$modification == 'UB') {
    data_f = data_f[Modifications %like% 'GlyGly']
  } else if (config$data$filters$modification == 'PH') {
    data_f = data_f[Modifications %like% 'Phospho']
  } else{
    cat(
      "\nERROR!!!! <<<",
      config$data$filters$modification,
      ">>> IS NOT A VALID config$data$filters$modification OPTION\n"
    )
    stop("CHECK HELP FOR VALID OPTIONS")
  }
  return(data_f)
}

# ------------------------------------------------------------------------------
#' @title Remove contaminants and empty proteins from the MaxQuant evidence file
#'
#' @description Remove contaminants and erronously identified 'reverse'
#' sequences by MaxQuant, in addition to empty protein ids
#' @param data (data.frame) of the Evidence file
#' @param verbose (logical) `TRUE` (default) shows function messages
#' @return (data.frame) without REV__ and CON__ Protein ids
#' @keywords cleanup, contaminants
#' @examples
#' ef <- artms_filterEvidenceContaminants(data = artms_data_ph_evidence)
#' @export
artms_filterEvidenceContaminants <- function(data,
                                             verbose = TRUE) {
  # Remove contaminants and reversed sequences (labeled by MaxQuant)
  data_selected <-
    data[grep("CON__|REV__", data$Proteins, invert = TRUE), ]
  # Remove empty proteins names
  blank.idx <- which(data_selected$Proteins == "")
  if (length(blank.idx) > 0)
    data_selected = data_selected[-blank.idx, ]
  if(verbose) cat(">> CONTAMINANTS CON__|REV__ REMOVED\n")
  return(data_selected)
}

# ------------------------------------------------------------------------------
#' @title Merge evidence.txt (or summary.txt) with keys.txt files 
#' @description Merge the evidence and keys files on the given columns
#' @param data (data.frame or char) The evidence data, either as data.frame or
#' the file name (and path). It also works for the summary.txt file
#' @param keys The keys data, either as a data.frame or file name (and path)
#' @param by (vector) specifying the columns use to merge the evidence and keys.
#' Default: `by=c('RawFile')`
#' @param isSummary (logical) TRUE or FALSE (default)
#' @param verbose (logical) `TRUE` (default) shows function messages
#' @return (data.frame) with the evidence and keys merged
#' @keywords merge, evidence, summary, keys
#' @examples
#' evidenceKeys <- artms_mergeEvidenceAndKeys(data = artms_data_ph_evidence,
#'                                            keys = artms_data_ph_keys)
#' @export
artms_mergeEvidenceAndKeys <- function(data, 
                                       keys, 
                                       by = c('RawFile'),
                                       isSummary = FALSE,
                                       verbose = TRUE) {

  if(verbose){
    cat(">> MERGING FILES\n")
    cat("\tIt might take a long time 
        (depending on the size of the evidence file)\n")
  }

  data <- .artms_checkIfFile(data)
  keys <- .artms_checkIfFile(keys)
  
  data <- .artms_checkRawFileColumnName(data)
  keys <- .artms_checkRawFileColumnName(keys)
  
  if(any(grepl("Experiment", colnames(keys)))){
    keys <- artms_changeColumnName(keys, "Experiment", "ExperimentKeys")
  }
  
  if(isSummary){
    if(any(grepl("Experiment", colnames(data)))){
      data <- subset(data, Experiment != "") 
    }
  }
  
  requiredColumns <- c('RawFile',
                   'IsotopeLabelType',
                   'Condition',
                   'BioReplicate', 
                   'Run')
  # Check that the keys file is correct
  if (any(!requiredColumns %in% colnames(keys))) {
    stop('Column names in keys not conform to schema. Required columns:\n', 
           sprintf('\t%s\n', requiredColumns))
  }
  
  # Check if the number of RawFiles is the same.
  unique_data <- sort(unique(data$RawFile))
  unique_keys <- sort(unique(keys$RawFile))
  
  if (length(unique_keys) != length(unique_data)) {
    keys_not_found <- setdiff(unique_keys, unique_data)
    data_not_found <- setdiff(unique_data, unique_keys)
    if(length(keys_not_found) != 0){
      cat(
        sprintf(
          "\tRaw data not in data file: %s\n",
          paste(keys_not_found, collapse = '\t'))
      )
    }
    if(length(data_not_found) != 0){
      cat(
        sprintf(
          "\tRaw data not in keys file: %s\n",
          paste(data_not_found, collapse = '\t')
        )
      )
    }
  }
  
  data <- merge(data, keys, by = by)
  return(data)
}


# ------------------------------------------------------------------------------
#' @title Convert the SILAC evidence file to MSstats format
#'
#' @description Converting the evidence file from a SILAC search to a format
#' compatible with MSstats. It basically modifies the Raw.files adding the
#' Heavy and Light label
#' @param evidence_file (char) Text filepath to the evidence file
#' @param output (char) Text filepath of the output name
#' @return (data.frame) with SILAC data processed for MSstats (and output file)
#' @keywords convert, silac, evidence
#' @examples \donttest{
#' evidence2silac <- artms_SILACtoLong(evidence_file = "silac.evicence.txt",
#'                                    output = "silac-evidence.txt")
#' }
#' @export
artms_SILACtoLong <- function(evidence_file, output) {
  file <- Sys.glob(evidence_file)
  cat(sprintf('>> PROCESSING SILAC EVIDENCE FILE\n'))
  tmp <- fread(file, integer64 = 'double')
  
  # reshape the data and split the heavy and light data
  tmp_long <- reshape2::melt(tmp, measure.vars = c('Intensity L', 
                                                   'Intensity H'))
  tmp_long[, Intensity := NULL]
  setnames(tmp_long, 'value', 'Intensity')
  setnames(tmp_long, 'variable', 'IsotopeLabelType')
  setnames(tmp_long, 'Raw file', 'Raw.file')
  levels(tmp_long$IsotopeLabelType) = c('L', 'H')
  tmp_long[!is.na(tmp_long$Intensity) &&
             tmp_long$Intensity < 1, ]$Intensity = NA
  write.table(
    tmp_long,
    file = output,
    sep = '\t',
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  cat("--- File ", output, " is ready\n")
  return(tmp_long)
}

# ------------------------------------------------------------------------------
# @title Pretty Labels for Heatmaps
#
# @description Generates pretty labels for the heatmaps.
# @param uniprot_acs (char) Uniprot accession id
# @param uniprot_ids (char) Uniprot entry id
# @param gene_names (car) Gene symbol
# @return Pretty labels for a heatmap
# @keywords internal, plots, pretty
.artms_prettyPrintHeatmapLabels <-
  function(uniprot_acs, uniprot_ids, gene_names) {
    result = paste(uniprot_acs, uniprot_ids, gene_names, sep = ' ')
    return(result)
  }

# ------------------------------------------------------------------------------
# @title Remove protein groups
#
# @description Remove the group of proteins ids separated by separated by `;`
# @param data (data.frame) with a `Proteins` column.
# @return (data.frame) with the protein groups removed
# @keywords maxquant, remove, proteingroups
.artms_removeMaxQProteinGroups <- function(data) {
  data_selected = data[grep(";", data$Proteins, invert = TRUE), ]
  return(data_selected)
}


# ------------------------------------------------------------------------------
#' @title Reshape the MSstats results file from long to wide format
#'
#' @description Converts the normal MSStats results.txt file into "wide" format
#' where each row represents a unique protein's results, and each column
#' represents the comparison made by MSStats. The fold change and p-value
#' of each comparison will be its own column.
#' @param results_msstats (char) Input file name and location
#' (MSstats `results.txt` file)
#' @param output_file (char) Output file name and location
#' (e.g. `results-wide.txt`). If `NULL` (default) returns an
#' R object (data.frame)
#' @param select_pvalues (char) Either
#' - `pvalue` or
#' - `adjpvalue` (default)
#' @param species (char) Specie name for annotation purposes.
#' Check `?artmsMapUniprot2Entrez` to find out more about the
#' supported species (e.g `species = "human"`)
#' @return (output file tab delimited) reshaped file with unique protein ids
#' and as many columns log2fc and adj.pvalues as comparisons available
#' @keywords msstats, results, wide, reshape
#' @examples
#' ph_results_wide <- artms_resultsWide(
#'                          results_msstats = artms_data_ph_msstats_results,
#'                          output_file = NULL,
#'                          species = "human")
#' @export
artms_resultsWide <- function(results_msstats,
                              output_file = NULL,
                              select_pvalues = "adjpvalue",
                              species) {
  cat(">> RESHAPING MSSTATS RESULTS TO wide FORMAT\n")
  results_msstats <- .artms_checkIfFile(results_msstats)
  
  if (select_pvalues == "adjpvalue") {
    input_l <-
      reshape2::melt(data <-
                       results_msstats[, c('Protein', 
                                           'Label', 
                                           'log2FC', 
                                           'adj.pvalue')], 
                     id.vars = c('Protein', 'Label'))
  } else if (select_pvalues == "pvalue") {
    input_l <-
      reshape2::melt(data = results_msstats[, c('Protein', 
                                                'Label', 
                                                'log2FC', 
                                                'pvalue')], 
                     id.vars = c('Protein', 'Label'))
  }
  
  ## then cast to get combinations of LFCV/PVAl and Label as columns
  input_w <-
    data.table::dcast(Protein ~ Label + variable,
                      data = input_l,
                      value.var = c('value'))
  suppressMessages(input_w <-
                     artms_annotationUniprot(input_w, "Protein", species))
  if (!is.null(output_file)) {
    write.table(
      input_w,
      file = output_file,
      eol = '\n',
      sep = '\t',
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE
    )
    cat("--- Results wide are out!\n")
  } else{
    return(input_w)
  }
}

# ------------------------------------------------------------------------------
# @title Correlation heatmaps of all the individual features
# @description Correlation heatmap using intensity values across all the
# conditions
# @param data_w (data.frame) resulting from the `.artms_castMaxQToWidePTM`
# function
# @param keys (data.frame) of the keys
# @param config (yaml.object) Configuration object (yaml loaded)
# @return (pdf) A correlation heatmap (suffix `-heatmap.pdf`)
# @keywords internal, heatmap, intensity, comparisons
.artms_sampleCorrelationHeatmap <- function (data_w, keys, config) {
  mat = log2(data_w[, 4:ncol(data_w), with = FALSE])
  mat[is.na(mat)] = 0
  mat_cor = cor(mat, method = 'pearson', use = 'everything')
  ## we want to make informarive row names so order by 
  ## RawFile because that's how data_w is ordered
  ordered_keys = keys[with(keys, order(RawFile)), ] 
  mat_names = paste(ordered_keys$Condition,
                    ordered_keys$BioReplicate,
                    ordered_keys$Run)
  colnames(mat_cor) = mat_names
  rownames(mat_cor) = mat_names
  colors_pos = colorRampPalette(RColorBrewer::brewer.pal("Blues", n = 5))(10)
  colors_neg = rev(colorRampPalette(RColorBrewer::brewer.pal("Reds", n =
                                                               5))(10))
  colors_tot = c(colors_neg, colors_pos)
  pheatmap(
    mat = mat_cor,
    cellwidth = 10,
    cellheight = 10,
    scale = 'none',
    filename = gsub('.txt', '-heatmap.pdf', config$files$output),
    color = colors_tot,
    breaks = seq(from = -1, to = 1, by = .1),
    fontfamily = "mono"
  )
  
}

# ------------------------------------------------------------------------------
# @title Barplot of peptide counts per biological replicate
#
# @description Total number of unique peptide identified per biological
# replicate
# @param data_f (char) Evidence file (same structure as the original)
# @param config (yaml.object) Configuration object
# @return (pdf) Barplot of peptide counts
# @keywords barplot, counts, peptides
.artms_samplePeptideBarplot <- function(data_f, config) {
  # set up data into ggplot compatible format
  data_f <-
    data.table(data_f,
               labels = paste(data_f$RawFile, 
                              data_f$Condition, 
                              data_f$BioReplicate))
  data_f <- data_f[with(data_f, order(labels, decreasing = TRUE)), ]
  
  # plot the peptide counts for all the samples TOGETHER
  p <- ggplot(data = data_f, aes(x = labels))
  p <-
    p + geom_bar() + theme(axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      family = 'mono'
    )) + ggtitle('Unique peptides per run\n after filtering') + coord_flip()
  ggsave(
    filename = gsub('.txt', '-peptidecounts.pdf', config$files$output),
    plot = p,
    width = 8,
    height = 10
  )
  
  w <- 10
  h <-
    ceiling((7 / 5 + 2) * ceiling(length(unique(data_f$Condition)) / 5))
  # plot the peptide counts for all the samples PER BAIT
  p <- ggplot(data = data_f, aes(x = as.factor(BioReplicate)))
  p <-
    p + geom_bar() + theme(axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      family = 'mono' )) + 
    ggtitle('Unique peptides per run\n after filtering') + 
    facet_wrap( ~ Condition, scales = 'free', ncol = 5)  + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(
    filename = gsub('.txt', '-peptidecounts-perBait.pdf', config$files$output),
    plot = p,
    width = w,
    height = h
  )
  
}

# ------------------------------------------------------------------------------
# @title Select significant hits
#
# @description Filtered data.frame with significant values (log2fc > 2 |
# log2fc < -2; adj.pvalue < 0.05) from the MSstats results
# @param mss_results (data.frame) of MSstats results
# @param labels (vector) of selected labels. Default: all (`*`)
# @param LFC (vector, int) with the negative and positive threshold. Default:
# c(-2, 2)
# @param whatPvalue (char) `pvalue` or `adj.pvalue` (default)?
# @param FDR (int) false discovery rate (adj.pvalue) threshold. Default: 0.05
# @return (data.frame) only with significant hits
# @keywords internal, significant, selections
.artms_significantHits <- function(mss_results,
                                   labels = '*',
                                   LFC = c(-2, 2),
                                   whatPvalue = 'adj.pvalue',
                                   FDR = 0.05) {
  selected_results <- mss_results[grep(labels, mss_results$Label),]
  
  if (whatPvalue == "adj.pvalue") {
    significant_proteins <-
      selected_results[(
        !is.na(selected_results$log2FC) &
          selected_results$adj.pvalue <= FDR &
          (
            selected_results$log2FC >= LFC[2] |
              selected_results$log2FC <= LFC[1]
          )
      ), 'Protein']
  } else if (whatPvalue == "pvalue") {
    significant_proteins <-
      selected_results[(
        !is.na(selected_results$log2FC) &
          selected_results$pvalue <= FDR &
          (
            selected_results$log2FC >= LFC[2] |
              selected_results$log2FC <= LFC[1]
          )
      ), 'Protein']
  } else{
    stop("The whatPvalue argument is wrong. Valid options: pvalue or adj.pvalue")
  }
  
  significant_results <-
    selected_results[selected_results$Protein %in% significant_proteins,]
  return(significant_results)
}


# ------------------------------------------------------------------------------
#' @title Outputs the spectral counts from the MaxQuant evidence file.
#'
#' @description Outputs the spectral counts from the MaxQuant evidence file.
#' @param evidence_file (char) Maxquant evidence file or data object
#' @param keys_file (char) Keys file with the experimental design or data object
#' @param output_file (char) Output file name (add `.txt` extension).
#' If `NULL` (default) it returns a data.frame object
#' @return A txt file with biological replicates, protein id, and spectral
#' count columns
#' @keywords spectral_counts, evidence
#' @examples
#' summary_spectral_counts <- artms_spectralCounts(
#'                                  evidence_file = artms_data_ph_evidence,
#'                                  keys_file = artms_data_ph_keys)
#' @export
artms_spectralCounts <- function(evidence_file,
                                 keys_file,
                                 output_file = NULL) {
  cat(">> EXTRACTING SPECTRAL COUNTS FROM THE EVIDENCE FILE\n")
  
  data <- .artms_checkIfFile(evidence_file)
  keys <- .artms_checkIfFile(keys_file)
  
  data <- .artms_checkRawFileColumnName(data)
  keys <- .artms_checkRawFileColumnName(keys)
  
  
  data <- artms_mergeEvidenceAndKeys(data, keys, by = c('RawFile'))
  data_sel <-
    data[, c('Proteins',
             'Condition',
             'BioReplicate',
             'Run',
             'MS.MS.count')]
  data_sel <-
    artms_changeColumnName(data_sel, 'MS.MS.count', 'spectral_counts')
  data_sel <-
    aggregate(
      spectral_counts ~ Proteins + Condition + BioReplicate + Run,
      data = data_sel,
      FUN = sum
    )
  data_sel <-
    data.frame(
      data_sel,
      AllCondition = paste(
        data_sel$Condition,
        data_sel$BioReplicate,
        data_sel$Run,
        sep = '_'
      )
    )
  
  if (!is.null(output_file)) {
    write.table(
      data_sel[, c('AllCondition', 'Proteins', 'spectral_counts')],
      file = output_file,
      eol = '\n',
      sep = '\t',
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE
    )
    cat(">> OUTPUT FILE <", output_file, "> is ready\n")
  } else{
    return(data_sel)
  }
}

# ------------------------------------------------------------------------------
# @title Remove white spaces
#
# @description Remove white spaces
# @param x (vector) A string
# @return (vector) with no white spaces
# @keywords internal, remove, whitespace
.artms_trim <- function (x) {
  gsub("^\\s+|\\s+$", "", x)
}

# ------------------------------------------------------------------------------
# @title Generate the contrast matrix required by MSstats from a txt file
# @description It simplifies the process of creating the contrast file
# @param contrast_file The text filepath of contrasts
# @param all_conditions a vector with all the conditions in the keys file
# @return (data.frame) with the contrast file in the format required by
# MSstats
# @author Tom Nguyen, David Jimenez-Morales
# @keywords check, contrast
.artms_writeContrast <- function(contrast_file, all_conditions = NULL) {
    input_contrasts <- readLines(contrast_file, warn = FALSE)
    #remove empty lines
    input_contrasts <-
      input_contrasts[vapply(input_contrasts, nchar, FUN.VALUE = 0) > 0]
    
    # check if contrast_file is old-format (i.e the contrast_file is a matrix)
    headers <- unlist(strsplit(input_contrasts[1], split = "\t"))
    if (length(headers) > 1) {
      newinput_contrasts <- c()
      for (i in 2:length(input_contrasts)) {
        newinput_contrasts <-
          c(newinput_contrasts, unlist(strsplit(input_contrasts[i], 
                                                split = "\t"))[1])
      }
      input_contrasts <- newinput_contrasts
    }
    
    # validate the input
    input_contrasts <- trimws(input_contrasts)
    valid <- TRUE
    accepted_chars <- c(LETTERS, letters, 0:9, '-', '_')
    for (x in input_contrasts) {
      if (x != "") {
        characs <- unlist(strsplit(x, split = ''))
        not_allowed_count <-
          length(which(!(characs %in% accepted_chars)))
        if (not_allowed_count > 0) {
          valid <- FALSE
          stop(paste(x, "is not a valid input"))
        }
        
        dash_count <- length(which(characs == '-'))
        if (dash_count != 1) {
          valid <- FALSE
          stop(paste(x, "needs to contain exactly 1 '-'"))
        }
      }
    }
    
    if (valid) {
      mat <- t(as.data.frame(strsplit(input_contrasts, split = '-')))
      rownames(mat) <- NULL
      conds <- sort(unique(c(mat[, 1], mat[, 2])))
      contrast_matrix <-
        matrix(0, nrow = nrow(mat), ncol = length(conds))
      colnames(contrast_matrix) <- conds
      rownames(contrast_matrix) <- input_contrasts
      
      for (i in seq_len(nrow(mat)) ) {
        cond1 <- mat[i, 1]
        cond2 <- mat[i, 2]
        contrast_matrix[i, cond1] <- 1
        contrast_matrix[i, cond2] <- -1
      }
      
      # check if conditions are all found in Evidence/Key
      if (!is.null(all_conditions)) {
        d <- setdiff(conds, all_conditions)
        if (length(d) > 0) {
          msg <-
            paste("These conditions are not found in the dataset:",
                  paste(d, collapse = ","))
          stop(msg)
        }
      }
      return (contrast_matrix)
    } else{
      stop(
        'Something went wrong while generating the contrast file.
        Please, let the developers know at <artms.help@gmail.com>'
      )
    }
}
