# ------------------------------------------------------------------------------
#' @title MaxQuant evidence file to SAINTexpress format
#'
#' @description Converts the MaxQuant evidence file to the 3 required files
#' by SAINTexpress. One can choose to either use the `spectral counts`
#' (use `msspc`) or the `intensities` (use `msint`) for the analysis.
#' @param evidence_file (char) The evidence file path and name
#' @param keys_file (char) Keys file with a SAINT column 
#' specifying test (`T`) and control (`C`) conditions
#' @param output_file (char) Output file name (must have extension .txt)
#' @param ref_proteome_file (char) Reference proteome path file name in
#' fasta format
#' @param quant_variable (char) choose either
#' - `msspc` (spectral counts, default) or
#' - `msint` (MS Intensity)
#' @param verbose (logical) `TRUE` (default) shows function messages
#' @return The 3 required files by SAINTexpress:
#' - `interactions.txt`
#' - `preys.txt`
#' - `baits.txt`
#' @keywords SAINT, SAINTexpress, APMS
#' @examples
#' # Testing that the files cannot be empty
#' artmsEvidenceToSaintExpress(evidence_file = NULL, 
#' keys_file = NULL, ref_proteome_file = NULL)
#' @export
artmsEvidenceToSaintExpress <- function(evidence_file,
                                           keys_file,
                                           ref_proteome_file,
                                           quant_variable = c('msspc','msint'),
                                           output_file, 
                                           verbose = TRUE) {
  if(verbose) message(">> CONVERTING TO SAINTexpress FORMAT ")
  
  if(is.null(evidence_file) & is.null(keys_file) & is.null(ref_proteome_file)){
    return("The evidence_file, keys_file, and ref_proteome_file 
           must not be empty")
  }

  if(any(missing(evidence_file) | 
         missing(keys_file) |
         missing(ref_proteome_file) | 
         missing(output_file)))
    stop("Missed (one or many) required argument(s)
         Please, check the help of this function to find out more")
    
  if (!grepl(".txt", output_file)) {
    stop(
      "Argument <output_file> must have the extension '.txt'"
    )
  }

  if(!file.exists(ref_proteome_file)){
    stop("The file ", ref_proteome_file, " does not exist ")
  }
  
  x <- fread(evidence_file, integer64 = 'double')
  keys <- fread(keys_file, integer64 = 'double')
  
  # Check Raw.file column
  x <- .artms_checkRawFileColumnName(x)
  keys <- .artms_checkRawFileColumnName(keys)
  
  # Check MS.MS
  x <- .artms_checkMSMSColumnName(x)
  
  
  if(verbose) message('>> VERIFYING DATA AND KEYS ')
  if (any(
    !c(
      'RawFile',
      'IsotopeLabelType',
      'Condition',
      'BioReplicate',
      'Run',
      'SAINT'
    ) %in% colnames(keys)
  )) {
    stop(
      'colnames in keys not conform to schema
      \tRawFile\tIsotopeLabelType\tCondition\tBioReplicate\tRun\tSAINT '
    )
  }
  
  saint_baits <- keys[, c('BioReplicate', 'Condition', 'SAINT'), with = FALSE]
  
  x <- artmsMergeEvidenceAndKeys(x, 
                                 keys, 
                                 by = c('RawFile'),
                                 verbose = verbose)
  
  data_f <- artmsFilterEvidenceContaminants(x = x, verbose = verbose)
  data_f <- .artms_removeMaxQProteinGroups(data_f)
  
  quant_variable <- match.arg(quant_variable)
  if(verbose) message(">> AGGREGATING ON ", quant_variable, " VALUES... ")
  ## aggregate over technical replicates if necessary
  if (quant_variable == 'msspc') {
    setnames(data_f, 'MS/MS Count', 'spectral_counts')
    data_f_agg <-
      aggregate(
        spectral_counts ~ BioReplicate+Condition+Proteins+Sequence+Charge,
        data = data_f,
        FUN = max
      )
    data_f_agg <-
      aggregate(
        spectral_counts ~ BioReplicate + Condition + Proteins,
        data = data_f_agg,
        FUN = sum
      )
  } else if (quant_variable == 'msint') {
    data_f_agg <-
      aggregate(
        Intensity ~ BioReplicate + Condition + Proteins + Sequence + Charge,
        data = data_f,
        FUN = max
      )
    data_f_agg <-
      aggregate(Intensity ~ BioReplicate + Condition + Proteins,
                data = data_f_agg,
                FUN = sum)
  } else{
    stop(" Wrong value for variable to quantify. 
         Please use 'msspc' or 'msint'")
  }
  
  ## IP name, bait name, prey name, and spectral counts or intensity values
  saint_interactions <- data_f_agg
  
  ref_proteome <- read.fasta(
    file = ref_proteome_file,
    seqtype = "AA",
    as.string = TRUE,
    set.attributes = TRUE,
    legacy.mode = TRUE,
    seqonly = FALSE,
    strip.desc = FALSE
  )
  p_lengths <- c()
  p_names <- c()
  for (e in ref_proteome) {
    p_lengths <- c(p_lengths, nchar(e[1]))
    p_names <- c(p_names, attr(e, 'name'))
  }
  ref_table <- data.table(names = p_names, lengths = p_lengths)
  ref_table[, uniprot_ac := gsub('([a-z,0-9,A-Z]+\\|{1})([A-Z,0-9,\\_]+)(\\|[A-Z,a-z,0-9,_]+)',
                                 '\\2',
                                 names)]
  ref_table[, uniprot_id := gsub('([a-z,0-9,A-Z]+\\|{1})([a-z,0-9,A-Z]+\\|{1})([A-Z,a-z,0-9,_]+)',
                                 '\\3',
                                 names)]
  
  unique_preys <- data.table(uniprot_ac = unique(data_f_agg$Proteins))
  saint_preys <- ref_table[, c('uniprot_ac', 'lengths', 'uniprot_id'), 
                          with = FALSE]
  saint_preys <- merge(unique_preys, saint_preys, by = 'uniprot_ac', 
                      all.x = TRUE)
  missing_lengths <- nrow(saint_preys[is.na(saint_preys$uniprot_id), ])
  saint_preys[is.na(saint_preys$uniprot_id), ]$uniprot_id = saint_preys[is.na(saint_preys$uniprot_id), ]$uniprot_ac
  if (missing_lengths > 0) {
    if(verbose)     
      message(
      sprintf(
        "--- WARNING! COMPUTING %s MISSING LENGTHS WITH THE MEDIAN LENGTH FROM THE DATASET ",
        missing_lengths))
    saint_preys[is.na(saint_preys$lengths), ]$lengths = median(saint_preys$lengths, na.rm = TRUE)
  }
  
  ## WRITE
  write.table(
    saint_baits,
    file = gsub('.txt', '-saint-baits.txt', output_file),
    eol = '\n',
    sep = '\t',
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  write.table(
    saint_preys,
    file = gsub('.txt', '-saint-preys.txt', output_file),
    eol = '\n',
    sep = '\t',
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  write.table(
    saint_interactions,
    file = gsub('.txt', '-saint-interactions.txt', output_file),
    eol = '\n',
    sep = '\t',
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  if(verbose){
    message(">> OUTPUT FILES: ")
    message("--- ", gsub('.txt', '-saint-baits.txt', output_file))
    message("--- ", gsub('.txt', '-saint-preys.txt', output_file))
    message("--- ", gsub('.txt', '-saint-interactions.txt', output_file))
  }
}
