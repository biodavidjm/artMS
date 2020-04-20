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
# @title Check the `MS/MS Count` column name on the evidence data.table
#
# @description Address case issue with the MS/MS Count column name
# @param (data.frame) keys or evidence files
# @return (data.frame) with the `MS/MS Count` column name
# @keywords internal msmscount, columname
.artms_checkMSMSColumnName <- function(df) {
  if (!('MS/MS Count' %in% colnames(df))) {
    if ("MS/MS count" %in% colnames(df)) {
      df <- artmsChangeColumnName(df, 'MS/MS count', 'MS/MS Count')
    } else{
      stop("cannot find the <MS/MS Count> column")
    }
  }
  return(df)
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
      df <- artmsChangeColumnName(df, 'Raw.file', 'RawFile')
    } else if ("Raw file" %in% colnames(df)) {
      df <- artmsChangeColumnName(df, 'Raw file', 'RawFile')
    } else{
      stop("cannot find the <Raw.file> column")
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
#' artms_data_ph_evidence <- artmsChangeColumnName(
#'                                dataset = artms_data_ph_evidence,
#'                                oldname = "Phospho..STY.",
#'                                newname = "PH_STY")
#' @export
artmsChangeColumnName <- function(dataset, oldname, newname) {
  if (!(oldname %in% colnames(dataset))) {
    stop(" The Column name provided <",
         oldname,
         "> was not found in the object provided ")
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
# @param x (data.frame) Evidence file
# @param config (yaml.object) Configuration object (opened yaml file)
# @param verbose (logical) `TRUE` (default) shows function messages
# @return (data.frame) filtered according to the options selected
# @keywords internal, filtering, remove, proteingroups, ptms
.artms_filterData <- function(x, 
                              config,
                              verbose = TRUE) {
  if(verbose) message(">> FILTERING ")
  
  if (config$data$filters$contaminants) {
    x <- artmsFilterEvidenceContaminants(x)
  }
  
  if (config$data$filters$protein_groups == 'remove') {
    if(verbose) message("-- READY TO REMOVE PROTEIN GROUPS: ")
    
    # SELECT FIRST THE LEADING RAZOR PROTEIN AS PROTEINS, DEPENDING ON THE 
    # MAXQUANT VERSION
    
    # Address the old version of maxquant
    if ( "Leading.Razor.Protein" %in% colnames(x) ) {
      x <- artmsChangeColumnName(x, "Leading.Razor.Protein", "Leading.razor.protein")
    }
    
    x <- artmsLeaveOnlyUniprotEntryID(x, "Proteins")
    x <- artmsLeaveOnlyUniprotEntryID(x, "Leading.razor.protein")
    
    # Check: if neither old version nor new version of leading razor protein
    # is found... stop it
    if ( "Leading.razor.protein" %in% colnames(x) ) {
      x$Proteins <- NULL
      data_f <- artmsChangeColumnName(x, "Leading.razor.protein", "Proteins")
      if(verbose) message("---- Use <Leading.razor.protein> as Protein ID")
    }else{
      stop("<Leading razor protein> column not found. Proteins groups cannot be removed")
    }
    
    # This is not necessary for now: the leading razor protein is unique
    # data_f <- .artms_removeMaxQProteinGroups(x)

  } else if (config$data$filters$protein_groups == 'keep') {
    if(verbose) message("-- PROTEIN GROUPS KEPT")
    data_f <- x
  } else{
    stop(
      "filtering option for <protein_groups> not valid 
      (options available: keep or remove)"
    )
  }

  
  # DEAL WITH OLD CONFIGURATION FILES WHEN config$data$filters$modification 
  # COULD BE EMPTY
  if (is.null(config$data$filters$modification)) {
    if(verbose) message("--- NO config$data$filters$modification provided. 
        Using 'AB' as default ")
  } else if (config$data$filters$modification == 'AB' |
             config$data$filters$modification == 'APMS') {
    if(verbose) message(sprintf("-- PROCESSING %s",
                                config$data$filters$modification))
  } else if (config$data$filters$modification == 'UB') {
    data_f = data_f[Modifications %like% 'GlyGly']
  } else if (config$data$filters$modification == 'PH') {
    data_f = data_f[Modifications %like% 'Phospho']
  } else if (config$data$filters$modification == 'AC') {
    data_f = data_f[Modifications %like% 'Acetyl']
  } else{
    stop("The config > data > filters > modification ",
         config$data$filters$modification," is not valid option")
  }
  return(data_f)
}

# ------------------------------------------------------------------------------
#' @title Remove contaminants and empty proteins from the MaxQuant evidence file
#'
#' @description Remove contaminants and erronously identified 'reverse'
#' sequences by MaxQuant, in addition to empty protein ids
#' @param x (data.frame) of the Evidence file
#' @param verbose (logical) `TRUE` (default) shows function messages
#' @return (data.frame) without REV__ and CON__ Protein ids
#' @keywords cleanup, contaminants
#' @examples
#' ef <- artmsFilterEvidenceContaminants(x = artms_data_ph_evidence)
#' @export
artmsFilterEvidenceContaminants <- function(x,
                                            verbose = TRUE) {
  # Remove contaminants and reversed sequences (labeled by MaxQuant)
  
  x <- .artms_checkIfFile(x)
  x <- .artms_checkRawFileColumnName(x)
  
  data_selected <- x[grep("CON__|REV__", x$Proteins, invert = TRUE), ]
  
  # Remove empty proteins names
  blank.idx <- which(data_selected$Proteins == "")
  if (length(blank.idx) > 0)
    data_selected = data_selected[-blank.idx, ]
  
  if(verbose) message("-- CONTAMINANTS CON__|REV__ REMOVED ")
  return(data_selected)
}

# ------------------------------------------------------------------------------
#' @title Merge evidence.txt (or summary.txt) with keys.txt files 
#' @description Merge the evidence and keys files on the given columns
#' @param x (data.frame or char) The evidence data, either as data.frame or
#' the file name (and path). It also works for the summary.txt file
#' @param keys The keys data, either as a data.frame or file name (and path)
#' @param by (vector) specifying the columns use to merge the evidence and keys.
#' Default: `by=c('RawFile')`
#' @param isSummary (logical) TRUE or FALSE (default)
#' @param verbose (logical) `TRUE` (default) shows function messages
#' @return (data.frame) with the evidence and keys merged
#' @keywords merge, evidence, summary, keys
#' @examples
#' evidenceKeys <- artmsMergeEvidenceAndKeys(x = artms_data_ph_evidence,
#'                                            keys = artms_data_ph_keys)
#' @export
artmsMergeEvidenceAndKeys <- function(x, 
                                     keys, 
                                     by = c('RawFile'),
                                     isSummary = FALSE,
                                     verbose = TRUE) {

  if(verbose){
    message(">> MERGING FILES ")
  }

  x <- .artms_checkIfFile(x)
  keys <- .artms_checkIfFile(keys)
  
  x <- .artms_checkRawFileColumnName(x)
  keys <- .artms_checkRawFileColumnName(keys)
  
  # Make sure that the Intensity column is not empty
  if(!isSummary){
    if(all(is.na(x$Intensity))){
      stop("The <Intensity> column of the evidence file is empty. artMS cannot continue")
    }  
  }
  
  if(any(grepl("Experiment", colnames(keys)))){
    keys <- artmsChangeColumnName(keys, "Experiment", "ExperimentKeys")
  }
  
  if(isSummary){
    if(any(grepl("Experiment", colnames(x)))){
      x <- subset(x, Experiment != "") 
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
           sprintf('\t%s ', requiredColumns))
  }
  
  # Check if the number of RawFiles is the same.
  unique_data <- sort(unique(x$RawFile))
  unique_keys <- sort(unique(keys$RawFile))
  
  if (length(unique_keys) != length(unique_data)) {
    keys_not_found <- setdiff(unique_keys, unique_data)
    data_not_found <- setdiff(unique_data, unique_keys)
    if(length(keys_not_found) != 0){
      message(
        sprintf(
          "--(-) Raw.files in keys not found in evidence file: %s\n",
          paste(keys_not_found, collapse = '\t'))
      ) 
    }
    if(length(data_not_found) != 0){
      if (!any(grepl("Total", data_not_found))){
        message(
          sprintf(
            "--(-) Raw.files in evidence not found in keys file: %s\n",
            paste(data_not_found, collapse = '\t')
          )
        )
      }
    }
  }
  
  x <- merge(x, keys, by = by)
  
  # Make the 0 values, NA values
  if(any(x$Intensity == 0, na.rm = TRUE)){
    zero_values <- length(x$Intensity[(x$Intensity == 0)])
    total_values <- length(x$Intensity)
    x$Intensity[(x$Intensity == 0)] <- NA
    message("---- (r) ", zero_values, " peptides with Intensity = 0 out of ", total_values, " replaced by NAs")
  }
  
  return(x)
}


# ------------------------------------------------------------------------------
#' @title Convert the SILAC evidence file to MSstats format
#'
#' @description Converting the evidence file from a SILAC search to a format
#' compatible with MSstats. It basically modifies the Raw.files adding the
#' Heavy and Light label
#' @param evidence_file (char) Text filepath to the evidence file
#' @param output (char) Text filepath of the output name. If NULL it does not
#' write the output
#' @param verbose (logical) `TRUE` (default) shows function messages
#' @return (data.frame) with SILAC data processed for MSstats (and output file)
#' @keywords convert, silac, evidence
#' @examples \donttest{
#' evidence2silac <- artmsSILACtoLong(evidence_file = "silac.evicence.txt",
#'                                    output = "silac-evidence.txt")
#' }
#' @export
artmsSILACtoLong <- function(evidence_file, 
                             output = NULL,
                             verbose = TRUE) {
  
  file <- Sys.glob(evidence_file)
  if(verbose) message(sprintf('>> PROCESSING SILAC EVIDENCE FILE '))
  tmp <- fread(file, integer64 = 'double')
  
  # reshape the data and split the heavy and light data
  tmp_long <- data.table::melt(tmp, 
                                measure.vars = c("Intensity L", "Intensity H"))
  
  tmp_long[, Intensity := NULL]
  setnames(tmp_long, 'value', 'Intensity')
  setnames(tmp_long, 'variable', 'IsotopeLabelType')
  setnames(tmp_long, 'Raw file', 'Raw.file')
  levels(tmp_long$IsotopeLabelType) = c('L', 'H')
  tmp_long[!is.na(tmp_long$Intensity) &&
             tmp_long$Intensity < 1, ]$Intensity = NA
  
  if(!is.null(output)){
    write.table(
      tmp_long,
      file = output,
      sep = '\t',
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE
    )
    if(verbose) message("--- File ", output, " is ready ")
  }
  colnames(tmp_long) <- gsub(" ", ".", colnames(tmp_long))
  colnames(tmp_long) <- gsub("/", ".", colnames(tmp_long))
  colnames(tmp_long) <- gsub("\\(", ".", colnames(tmp_long))
  colnames(tmp_long) <- gsub("\\)", ".", colnames(tmp_long))
  tmp_long <- as.data.frame(tmp_long)
  return(tmp_long)
}


# ------------------------------------------------------------------------------
# @title Merge keys and Evidence from SILAC experiments
#
# @description Merge keys and Evidence from SILAC experiments
# @param evisilac (char) Output from artmsSILACtoLong
# @param keysilac (char) keys files with SILAC details
# @return df with both evidence and keys from silac merge
# @keywords internal, silac, merge
.artmsMergeSilacEvidenceKeys <- function(evisilac, 
                                         keysilac){
  
  
  evisilac <- .artms_checkIfFile(evisilac)
  evisilac <- .artms_checkRawFileColumnName(evisilac)
  
  keys <- .artms_checkIfFile(keysilac)
  keys <- .artms_checkRawFileColumnName(keys)
  
  # Check the labels from the keys file
  hlvalues <- unique(keys$IsotopeLabelType)
  hl2find <- c("L","H")
  
  if(length(setdiff(hlvalues, hl2find)) > 0){
    stop("The IsotopeLabelType available in the keys file are not from a 
         SILAC experiment: they must be H and L")
  }
  
  evisilac$RawFile = paste(evisilac$RawFile, 
                           evisilac$IsotopeLabelType, 
                           sep = '')
  keysilac$RawFile = paste(keysilac$RawFile, 
                       keysilac$IsotopeLabelType, 
                       sep = '')
  keysilac$Run = paste(keysilac$IsotopeLabelType, 
                   keysilac$Run , 
                   sep = '')
  
  evisilac$IsotopeLabelType = 'L'
  keysilac$IsotopeLabelType = 'L'
  
  df <- artmsMergeEvidenceAndKeys(evisilac, 
                                 keysilac, 
                                 by = c('RawFile', 'IsotopeLabelType'),
                                 verbose = FALSE)
  return(df)
}
  

# ------------------------------------------------------------------------------
# @title Remove protein groups
#
# @description Remove the group of proteins ids separated by separated by `;`
# @param x (data.frame) with a `Proteins` column.
# @return (data.frame) with the protein groups removed
# @keywords maxquant, remove, proteingroups
.artms_removeMaxQProteinGroups <- function(x) {
  data_selected = x[grep(";", x$Proteins, invert = TRUE), ]
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
#' @param verbose (logical) `TRUE` (default) shows function messages
#' @return (output file tab delimited) reshaped file with unique protein ids
#' and as many columns log2fc and adj.pvalues as comparisons available
#' @keywords msstats, results, wide, reshape
#' @examples
#' ph_results_wide <- artmsResultsWide(
#'                          results_msstats = artms_data_ph_msstats_results,
#'                          output_file = NULL,
#'                          species = "human")
#' @export
artmsResultsWide <- function(results_msstats,
                             output_file = NULL,
                             select_pvalues = c("adjpvalue", "pvalue"),
                             species,
                             verbose = TRUE) {
  
  if(any(missing(results_msstats) | 
         missing(species)))
    stop("Missed (one or many) required argument(s)
         Please, check the help of this function to find out more")
  
  if(verbose) message(">> RESHAPING MSSTATS RESULTS TO wide FORMAT ")
  results_msstats <- .artms_checkIfFile(results_msstats)
  
  select_pvalues <- match.arg(select_pvalues)
  pvals <- if(select_pvalues == "adjpvalue") "adj.pvalue" else "pvalue"
  selectedColumns <- c('Protein', 'Label', 'log2FC', pvals)
  input_l <- data.table::melt(data = results_msstats[,selectedColumns], 
                              id.vars = c('Protein', 'Label'))
  
  ## then cast to get combinations of LFCV/PVAl and Label as columns
  input_w <- data.table::dcast(Protein ~ Label + variable,
                               data = input_l,
                               value.var = c('value'))
  suppressMessages(input_w <- artmsAnnotationUniprot(input_w, 
                                                      "Protein", 
                                                      species))
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
    if(verbose) message("--- Results wide are out! ")
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
#' @title Outputs the spectral counts from the MaxQuant evidence file.
#'
#' @description Outputs the spectral counts from the MaxQuant evidence file.
#' @param evidence_file (char) Maxquant evidence file or data object
#' @param keys_file (char) Keys file with the experimental design or data object
#' @param output_file (char) Output file name (add `.txt` extension).
#' If `NULL` (default) it returns a data.frame object
#' @param verbose (logical) `TRUE` (default) shows function messages
#' @return A txt file with biological replicates, protein id, and spectral
#' count columns
#' @keywords spectral_counts, evidence
#' @examples
#' summary_spectral_counts <- artmsSpectralCounts(
#'                                  evidence_file = artms_data_ph_evidence,
#'                                  keys_file = artms_data_ph_keys)
#' @export
artmsSpectralCounts <- function(evidence_file,
                                 keys_file,
                                 output_file = NULL,
                                 verbose = TRUE) {
  if(verbose) message(">> EXTRACTING SPECTRAL COUNTS FROM THE EVIDENCE FILE ")
  
  x <- .artms_checkIfFile(evidence_file)
  keys <- .artms_checkIfFile(keys_file)
  
  x <- .artms_checkRawFileColumnName(x)
  keys <- .artms_checkRawFileColumnName(keys)
  
  
  x <- artmsMergeEvidenceAndKeys(x, 
                                 keys, 
                                 by = c('RawFile'),
                                 verbose = verbose)
  data_sel <-
    x[, c('Proteins',
             'Condition',
             'BioReplicate',
             'Run',
             'MS.MS.count')]
  data_sel <-
    artmsChangeColumnName(data_sel, 'MS.MS.count', 'spectral_counts')
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
    if(verbose) message(">> OUTPUT FILE <", output_file, "> is ready ")
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
.artms_writeContrast <- function(contrast_file, 
                                 all_conditions = NULL) {
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
          stop(paste(x, " is not a valid input"))
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
