# ------------------------------------------------------------------------------
#' @title MaxQuant evidence file to SAINTq format
#'
#' @description Converts the MaxQuant evidence file to the required files
#' by SAINTq. The user can choose to use either peptides with `spectral counts`
#' (use `msspc`) or the all the peptides (use `all`) for the analysis. 
#' The quantitative can be also chosen (either MS Intensity or Spectral Counts)
#' @param evidence_file (char or data.frame) The evidence file path and name, 
#' or data.frame
#' @param keys_file (char) Keys file with a SAINT column specifying
#' test (`T`) and control (`C`) conditions
#' @param output_dir (char) New directory to create and save files. 
#' Default is current directory (recommended to provide a new folder name).
#' @param sc_option (char). Filter peptides with spectral counts only.
#' Two options:
#' - `msspc`: use only peptides with spectral_counts
#' - `all` (default): all peptides detected (including the one resulting from 
#' the MaxQuant 'Match between run' algorithm)
#' @param quant_variable (char) Select the quantitative variable. 
#' Two options available:
#' - `msint`: MS Intensity (default)
#' - `msspc`: MS.MS.count (Spectral Counts)
#' @param fractions (logical) `TRUE` for 2D proteomics (fractions). 
#' Default: `FALSE`
#' @param verbose (logical) `TRUE` (default) shows function messages
#' @return The input files requires to run SAINTq
#' @details After running the script, the new specified folder should contain
#' the folling files:
#' - saintq-config-peptides
#' - saintq-config-proteins
#' - saintq_input_peptides.txt
#' - saintq_input_proteins.txt
#' 
#' Then `cd` into the new folder and run either of the following two options
#' (assuming that `saintq` is installed in your linux/unix/mac os x system):
#' 
#' `> saintq config-saintq-peptides`
#' 
#' or
#' 
#' `> saintq config-saintq-proteins`
#' @keywords SAINT, SAINTq, APMS
#' @examples
#' # Testing that the files cannot be empty
#' artmsEvidenceToSAINTq   (evidence_file = NULL, 
#'                                    keys_file = NULL, 
#'                                    output_dir = NULL)
#' @export
artmsEvidenceToSAINTq    <- function(evidence_file,
                                         keys_file,
                                         output_dir = ".",
                                         sc_option = c("all", "msspc"),
                                         fractions = FALSE,
                                         quant_variable = c('msint','msspc'),
                                         verbose = TRUE){

  if(verbose){
    message(">> GENERATING A SAINTq INPUT FILE ")
    message(">> CHECKING THE keys FILE FIRST ")
  }  
  
  if(is.null(evidence_file) & is.null(keys_file) & is.null(output_dir)){
    return("The evidence_file, keys_file and output_dir must not be NULL")
  }
  
  if(any(missing(evidence_file) | missing(keys_file)))
    stop("Missed (one or many) required argument(s)
         Please, check the help of this function to find out more")
  
  keys <- .artms_checkIfFile(keys_file)
  keys <- .artms_checkRawFileColumnName(keys)
  
  if(fractions){
    if(verbose) message("--- VERIFYING THAT THE INFORMATION ABOUT fractions IS AVAILABLE ")
    requiredColumns <- c('RawFile','IsotopeLabelType','Condition',
                         'BioReplicate','Run', 
                         'FractionKey', 'SAINT')
    if(any(! requiredColumns %in% colnames(keys)))
      stop('Column names in keys not conform to schema. Required columns:', 
           sprintf('\t%s ',requiredColumns))
  }else{
    requiredColumns <- c('RawFile','IsotopeLabelType','Condition',
                         'BioReplicate','Run', 'SAINT')
    if(any(! requiredColumns %in% colnames(keys)))
      stop('Column names in keys not conform to schema. Required columns:', 
           sprintf('\t%s ', requiredColumns))
  }
  
  # EVIDENCE:
  datamerged <- artmsMergeEvidenceAndKeys(evidence_file, 
                                           keys_file,
                                           verbose = verbose)
  
  # SELECTING THE Leading.razor.protein
  
  datamerged <- subset(datamerged, select = -Proteins)
  if( ('Leading.razor.protein' %in% colnames(datamerged)) ) {
    if(verbose) message('--- Making the <Leading.Razor.Protein> the <Proteins> column ')
    names(datamerged)[grep('Leading.razor.protein', 
                           names(datamerged))] <- 'Proteins'
  } else if('Leading.Razor.Protein' %in% colnames(datamerged) ) {
    if(verbose) message('--- Making the <Leading.Razor.Protein> the <Proteins> column ')
    names(datamerged)[grep('Leading.Razor.Protein', names(datamerged))] <-
      'Proteins'
  } else{
    stop("there is no <leading.razor.protein> column in this evidence file.")
  }
  
  datamerged$Proteins <- gsub("(sp\\|)(.*)(\\|.*)", "\\2", datamerged$Proteins )

  ## ONLY VALUES WITH SPECTRAL COUNTS
  sc_option <- match.arg(sc_option)
  if(sc_option == "msspc"){
    if(verbose) message("--- Selecting peptides with spectral count only ")
    before <- dim(datamerged)[1]
    datamerged <- datamerged[which(datamerged$MS.MS.Count > 0),]
    after <- dim(datamerged)[1]
    keepingPercent <- (after*100)/before
    if(verbose){
      message("\t+--> Before:", before," ")
      message("\t+--> After:", after," (Keeping:", keepingPercent,"%) ")
    }
  }else if(sc_option == "all"){
    if(verbose)
      message("--- ALL peptides with intensities will be used to generate the 
      saintq input file (indepependently of the number of spectral counts ")
  }else{stop("<sc_option> argument must be either 'sc' or 'all'")
  }

  # Remove empty proteins
  if(verbose) message("--- Removing empty protein ids (if any) ")
  if(length(which(datamerged$Proteins==""))>0){
    datamerged <- datamerged[-which(datamerged$Proteins==""),]
  }
  
  # Remove protein groups
  if(verbose) message("--- Removing Protein Groups (if any) ")
  datamerged <- .artms_removeMaxQProteinGroups(datamerged)
  
  # Removing Contaminants
  if(verbose) message("--- Removing contaminants")
  data_f2 <- artmsFilterEvidenceContaminants(datamerged,
                                              verbose = verbose)
  
  quant_variable <- match.arg(quant_variable)
  if(quant_variable == "msint"){
    # Set the intensity as numeric to avoid overflow problems
    data_f2$Intensity = as.numeric(data_f2$Intensity)
    data_f2 <- data_f2[!is.na(data_f2$Intensity),]
  }else if(quant_variable == "msspc"){
    # hack to use spectral counts instead of intensities
    data_f2$Intensity <- NA
    data_f2$Intensity <- data_f2$MS.MS.Count
    data_f2$Intensity = as.numeric(data_f2$Intensity)
    data_f2 <- data_f2[!is.na(data_f2$Intensity),]
  }
  
  # create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Combine all the fractions if this is a fractioning experiment by 
  # summing them up
  if (fractions){
    # Sum up all the fractions first
    data_f2_fa <- aggregate(
      Intensity~Sequence+Proteins+Condition+BioReplicate+Run, 
      data=data_f2, 
      FUN = sum)
  }else{
    data_f2_fa <- data_f2
  }

  # Take only unique values of both sequences and proteins.
  data_f2uniques <- unique(data_f2_fa[,c("Proteins", "Sequence")])
  protPep <- aggregate(Sequence ~ Proteins, 
                       data_f2uniques, 
                       FUN = paste, collapse="|" )
  
  # Using Intensities for both Peptides and Proteins
  protBiorepIntensity  <- data.table::dcast(
    data=data_f2_fa[,c("Proteins","BioReplicate","Intensity")], 
    Proteins~BioReplicate, 
    value.var = "Intensity", sum, na.rm = TRUE, fill=0 )
  peptideBiorepIntensity <- data.table::dcast(
    data=data_f2_fa[,c("Proteins","Sequence","BioReplicate","Intensity")], 
    Proteins+Sequence~BioReplicate, 
    value.var = "Intensity", fun.aggregate = sum, na.rm = TRUE, fill=0 )

  # PROTEINS: extra step to add the information about peptides
  almost <- merge(protBiorepIntensity, protPep, by="Proteins")
  final_result <- almost[,c(1,dim(almost)[2],2:(dim(almost)[2]-1)) ]
  
  # SAINTQ HEADER (two extra rows on top)
  rekey <- keys[keys$BioReplicate %in% data_f2_fa$BioReplicate,]
  
  if(fractions){
    rekey <- rekey[c('Condition', 'BioReplicate', 'SAINT')]
    rekey <- unique(rekey)
  }
  
  x <- t(rekey[,c('SAINT','Condition','BioReplicate')])
  extra <- cbind(c('','','Proteins'),c('','','Sequence'))
  header <- t(cbind(extra, x))
  theader <- t(header)
  checkthis <- data.frame(theader, row.names = NULL, stringsAsFactors = FALSE)
  names(checkthis) = checkthis[3,]
  
  # PROTEINS: ADDING HEADER, merging based on row names
  proteinssaintqheader <- rbind(checkthis, final_result)
  # SEQUENCES: Adding HEADER:
  sequencesaintqheader <- rbind(checkthis, peptideBiorepIntensity)
  
  # Writing the files for PROTEIN
  output <- paste0(output_dir,'/saintq_input_proteins.txt')
  write.table(proteinssaintqheader, 
              output, sep='\t',
              row.names = FALSE, col.names= FALSE, quote = FALSE)
  
  outconfig_protein <- paste0(output_dir,'/config-saintq-proteins')
  cat ("
       ### SAINTq parameter file
       ## use # to mark a line as comment
       
       ## normalize control intensities
       normalize_control=false
       
       ## name of file with intensities
       input_filename=saintq_input_proteins.txt
       
       ## valid: protein, peptide, fragment
       input_level=protein
       
       ## column names
       protein_colname=Proteins
       
       ## control bait selection rules
       compress_n_ctrl=100
       
       ## test bait replicate selection rules
       compress_n_rep=100
       ", file=outconfig_protein)
  
  
  # WRITING the files for SEQUENCES
  outsequences <- paste0(output_dir,'/saintq_input_peptides.txt')
  write.table(sequencesaintqheader, 
              outsequences, 
              sep = '\t', 
              row.names = FALSE, 
              col.names = FALSE, 
              quote = FALSE)
  
  outconfig_peptide <- paste0(output_dir,'/config-saintq-peptides')
  cat ("
       ### SAINTq parameter file
       ## use # to mark a line as comment
       
       ## normalize control intensities
       normalize_control=false
       
       ## name of file with intensities
       input_filename=saintq_input_peptides.txt
       
       ## type of intensity
       ## valid: protein, peptide, fragment
       input_level=peptide
       
       ## column names
       protein_colname=Proteins
       pep_colname=Sequence
       
       ## control bait selection rules
       compress_n_ctrl=100
       
       ## test bait replicate selection rules
       compress_n_rep=100
       
       ## peptide selection rules
       min_n_pep=3
       best_prop_pep=0.5
       
       ", file=outconfig_peptide)

  if(verbose){
    message(">> NEW 4 FILES CREATED: ")
    message("\t- saintq-config-peptides ")
    message("\t- saintq-config-proteins ")
    message("\t- saintq_input_peptides.txt ")
    message("\t- saintq_input_proteins.txt  ")
    message("--- Now get into the folder and run either:\n
        > saintq config-saintq-peptides\n
        or\n
        > saintq config-saintq-proteins")
    message(">> DONE! ")
  }
}
