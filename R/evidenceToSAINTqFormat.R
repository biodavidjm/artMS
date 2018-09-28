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
#' @param output_dir (char) New directory to create ans save files
#' @param sc_option (char). Filter peptides with spectral counts only.
#' Two options:
#' - `msspc`: use only peptides with spectral_counts
#' - `all` (default): all peptides detected (including the one resulting from 
#' the MaxQuant 'Match between run' algorithm)
#' @param quant_variable (char) Select the quantitative variable. 
#' Two options available:
#' - `msint`: MS Intensity
#' - `msspc`: MS.MS.count (Spectral Counts)
#' @param fractions (logical) `TRUE` for 2D proteomics (fractions). 
#' Default: `FALSE`
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
#' artms_evidenceToSAINTqFormat(evidence_file = NULL, 
#'                                    keys_file = NULL, 
#'                                    output_dir = NULL)
#' @export
artms_evidenceToSAINTqFormat <- function(evidence_file,
                                         keys_file,
                                         output_dir,
                                         sc_option = "all",
                                         fractions = FALSE,
                                         quant_variable = "msint"){
  
  cat(">> GENERATING A SAINTq INPUT FILE\n")
  
  cat(">> CHECKING THE keys FILE FIRST\n")
  
  if(is.null(evidence_file) & is.null(keys_file) & is.null(output_dir)){
    return("The evidence_file, keys_file and output_dir must not be empty")
  }
    
  keys <- .artms_checkIfFile(keys_file)
  keys <- .artms_checkRawFileColumnName(keys)
  
  if(fractions){
    cat("--- VERIFYING THAT THE INFORMATION ABOUT fractions IS AVAILABLE\n")
    if(any(!c('RawFile','IsotopeLabelType','Condition',
              'BioReplicate','Run', 
              'FractionKey', 'SAINT') %in% colnames(keys))){
      cat('\nERROR: COLUMN NAMES IN KEYS NOT CONFORM TO SCHEMA. 
One of these is lost\n
          \tRawFile\n\tIsotopeLabelType\n\tCondition\n\tBioReplicate
          \tRun\n\n\tFractionKey\n\n\tSAINT\n\n')
      stop('Please, try again once revised\n\n')
    }
  }else{
    if(any(!c('RawFile','IsotopeLabelType',
              'Condition','BioReplicate',
              'Run','SAINT') %in% colnames(keys))){
      cat('\nERROR: COLUMN NAMES IN KEYS NOT CONFORM TO SCHEMA. 
One of these is lost\n
\tRawFile\n\tIsotopeLabelType\n\tCondition\n\tBioReplicate\n\tRun\n\n\tSAINT\n')
      stop('Please, try again once revised\n\n')
    }
  }
  
  # EVIDENCE:
  datamerged <- artms_mergeEvidenceAndKeys(evidence_file, keys_file)
  
  # SELECTING THE Leading.razor.protein
  
  datamerged <- subset(datamerged, select = -Proteins)
  if( ('Leading.razor.protein' %in% colnames(datamerged)) ) {
    cat('--- Making the <Leading.Razor.Protein> the <Proteins> column\n')
    names(datamerged)[grep('Leading.razor.protein', names(datamerged))] <-
      'Proteins'
  } else if('Leading.Razor.Protein' %in% colnames(datamerged) ) {
    cat('--- Making the <Leading.Razor.Protein> the <Proteins> column\n')
    names(datamerged)[grep('Leading.Razor.Protein', names(datamerged))] <-
      'Proteins'
  } else{
    stop("\n\n\n\tOH NO! THERE IS NO Leading.razor.protein COLUMN IN THIS 
         EVIDENCE FILE!!\n\n\n")
  }
  
  datamerged$Proteins <- gsub("(sp\\|)(.*)(\\|.*)", "\\2", datamerged$Proteins )

  ## ONLY VALUES WITH SPECTRAL COUNTS
  
  if(sc_option == "msspc"){
    cat("--- Selecting peptides with spectral count only\n")
    before <- dim(datamerged)[1]
    datamerged <- datamerged[which(datamerged$MS.MS.Count > 0),]
    after <- dim(datamerged)[1]
    keepingPercent <- (after*100)/before
    cat("\t+--> Before:", before,"\n")
    cat("\t+--> After:", after," (Keeping:", keepingPercent,"%)\n")
  }else if(sc_option == "all"){
    cat("--- ALL peptides with intensities will be used to generate the 
        saintq input file (indepependently of the number of spectral counts\n")
  }else{
    cat("\n\nWAIT A MINUTE: sc_option MUST BE EITHER 'sc' or 'all'\n\n")
    stop("\n\nTRY AGAIN WHEN READY\n\n")
  }

  # Remove empty proteins
  cat("--- Removing empty protein ids (if any)\n")
  if(length(which(datamerged$Proteins==""))>0){
    datamerged <- datamerged[-which(datamerged$Proteins==""),]
  }
  
  # Remove protein groups
  cat("--- Removing Protein Groups (if any)\n")
  datamerged <- .artms_removeMaxQProteinGroups(datamerged)
  
  # Removing Contaminants
  cat("--- Removing contaminants")
  data_f2 <- artms_filterEvidenceContaminants(datamerged)
  
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
  
  # Combine sequence and charge
  # data_f2$Sequence <- paste0(data_f2$Sequence,"_",data_f2$Charge)
  
  # Create the directory
  dir.create(output_dir, showWarnings = FALSE)
  
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

  cat(">> FOLDER <",output_dir,"> SHOULD CONTAIN 4 FILES:\n")
  cat("\t- saintq-config-peptides\n")
  cat("\t- saintq-config-proteins\n")
  cat("\t- saintq_input_peptides.txt\n")
  cat("\t- saintq_input_proteins.txt\n\n")
  cat("--- Now get into the folder and run either:
> saintq config-saintq-peptides
or
> saintq config-saintq-proteins\n")
  cat(">> DONE!\n")
}
