# ------------------------------------------------------------------------------
#' @title Convert MaxQuant Evidence file into a Protein Prospector like format 
#' to run through the MIST pipeline using spectral counts.
#' 
#' @description Converts MaxQuant evidence file into a file format compatible 
#' with the MiST pipeline using `MS/MS.Count`. Note that this is the 
#' MiST *data* file, and that an additional *keys* file will have to be 
#' constructed before running MiST. Multiple species can be searched at once, 
#' simply separate them by a "-". (eg. `HUMAN-MOUSE`)
#' @param input_file (char) MaxQuant evidence file and location
#' @param keys_file (char) Keys file with the experimental details
#' @param quant_variable (char) Select the quantitative variable. Two options available:
#' - `msint`: MS Intensity
#' - `msspc`: MS.MS.count (Spectral Counts)
#' @param output_file (char) Output file name
#' @param species (char) Specie name. If several, used a `dash` symbol to separate them
#' @param uniprot_dir Directory with the uniprot files with the mapping 
#' information. Default '~/Box Sync/db/mist/'
#' @return MIST compatible input files (reformatted evidence and keys)
#' @keywords mist, evidence, keys, apms
#' @examples \donttest{
#' artms_evidenceToMISTformat(input_file = "a549-PB1-evidence.txt", 
#'                            keys_file = "a549-PB1-keys.txt", 
#'                            quant_variable = "msint", 
#'                            output_file = "a549-PB1-mist.txt",
#'                            species = "HUMAN-FLUOMICS",
#'                            uniprot_dir = "~/Box Sync/db/mist/")
#' }
#' @export
artms_evidenceToMISTformat <- function(input_file, quant_variable, keys_file, output_file, species, uniprot_dir){
  cat('\n>> GENERATING INPUT FILES FOR MIST\n')
  if(quant_variable == "msint"){
    cat("---USING MS Intensity\n")
  }else if(quant_variable == "msspc"){
    cat("---USING MS.MS.Count (spectral counts)\n")
  }else{
    stop("\n << ",quant_variable," >> NOT ALLOWED. ONLY msint OR msspc ALLOWED\n")
  }
  cat('\tREADING IN DATA AND KEYS\n')
  
  if(is.null(species)){
    species <- "HUMAN"
    cat("\t---species NOT SPECIFIED. USING DEFAULT: ",species,"\n")
  }
  
  if(is.null(uniprot_dir)){
    uniprot_dir <- '~/Box Sync/db/mist/'
    cat("\t---uniprot_dir NOT SPECIFIED. USING DEFAULT:  ",uniprot_dir,"\n")
  }
  
  data <- data.table(read.delim(input_file, stringsAsFactors= FALSE))
  keys <- data.table(read.delim(keys_file, stringsAsFactors = FALSE))
  
  cat('\tCHECKING DATA AND KEYS COLUMN NAMES\n')
  tryCatch(setnames(data, 'Raw file', 'RawFile'), error=function(e) cat('\t---Raw file in evidence not found: trying Raw.file instead\n'))
  tryCatch(setnames(data, 'Raw.file', 'RawFile'), error=function(e) cat('\t---Raw.file in evidence not found: trying Raw file instead\n'))
  tryCatch(setnames(keys, 'Raw file', 'RawFile'), error=function(e) cat('\t---Raw file in keys not found: trying Raw.file instead\n'))
  tryCatch(setnames(keys, 'Raw.file', 'RawFile'), error=function(e) cat('\t---Raw.file in keys not found: trying Raw file instead\n'))
  
  if(quant_variable == "msspc"){
    tryCatch(setnames(data,'MS/MS Count','ms_spectral_counts'), error=function(e) cat('\t---MS/MS Count column not found in evidence file: trying MS.MS.Count instead\n'))
    tryCatch(setnames(data,'MS.MS.Count','ms_spectral_counts'), error=function(e) cat('\t---MS.MS.Count column not found in evidence file: trying MS/MS Count instead\n'))
    tryCatch(setnames(data,'MS.MS.count','ms_spectral_counts'), error=function(e) cat('\t---MS.MS.count column not found in evidence file: trying MS/MS Count instead\n'))
  }else if(quant_variable == "msint"){
    tryCatch(setnames(data,'Intensity','ms_intensity'), error=function(e) stop('\t--- INTENSITY NOT FOUND IN THE evidence FILE!!\n\n'))
  }else{
    stop("\n << ",quant_variable," >> NOT ALLOWED. ONLY msint OR msspc ALLOWED\n")
  }
  
  
  cat('\n\tVERIFYING DATA AND KEYS\n')
  if(!'IsotopeLabelType' %in% colnames(data)) data[,IsotopeLabelType:='L']
  data <- artms_mergeEvidenceAndKeys(data, keys, by = c('RawFile','IsotopeLabelType'))

  if(quant_variable == "msspc"){
    data_sel <- data[,c('Proteins','Condition','BioReplicate','Run','RawFile','ms_spectral_counts'),with= FALSE]
    data_sel <- aggregate( ms_spectral_counts ~ Proteins+Condition+BioReplicate+Run+RawFile, data=data_sel, FUN = sum)
  }else if(quant_variable == "msint"){
    data_sel <- data[,c('Proteins','Condition','BioReplicate','Run','RawFile','ms_intensity'),with= FALSE]
    data_sel <- aggregate( ms_intensity ~ Proteins+Condition+BioReplicate+Run+RawFile, data=data_sel, FUN = sum)
  }else{
    stop("\n << ",quant_variable," >> NOT ALLOWED. ONLY msint OR msspc ALLOWED\n")
  }
  
  data_sel <- data.frame(data_sel, bait_name=paste(data_sel$Condition, data_sel$BioReplicate, data_sel$Run, sep='_'))
  
  # clean up proteins & annotate
  #~~~~~~~~~~~~~~~~~~~~~~~
  # remove CON's
  if( length(grep("^CON__",data_sel$Proteins))>0 ) data_sel = data_sel[-grep("^CON__",data_sel$Proteins),]
  if( length(grep("^REV__",data_sel$Proteins))>0 ) data_sel = data_sel[-grep("^REV__",data_sel$Proteins),]
  # remove the party sets
  if( length(grep(";",data_sel$Proteins))>0 ) data_sel = data_sel[-grep(";",data_sel$Proteins),]     # NOTE!!! We lose a lot of entries this way... :\
  # keep only uniprot id
  #data_sel$uniprot_id = gsub("^.*\\|","", data_sel$Proteins)
  data_sel$Proteins = gsub("(^.*\\|)([A-Z0-9]+)(\\|.*$)","\\2",data_sel$Proteins)
  # remove blank protein names
  if(any(data_sel$Proteins == "")){ data_sel <- data_sel[-which(data_sel$Proteins == ""),]}
  
  # Add this column that it is required for the preprocessing done by mist 
  # (this column is generated by Prospector and it is used to eleminate peptides. In this case, does not apply but it has to be there)
  data_sel$ms_unique_pep = ""
  # re-order

  if(quant_variable == "msspc"){
    data_sel <- data_sel[,c("RawFile",'Proteins','ms_unique_pep', 'ms_spectral_counts')]
    # RENAMING!
    names(data_sel) = c('id','ms_uniprot_ac','ms_unique_pep','ms_spectral_counts')
    # remove interactions with ms_spectral_counts=0
    if(any(data_sel$ms_spectral_counts == 0)) { data_sel <- data_sel[-which(data_sel$ms_spectral_counts==0),]}
  }else if(quant_variable == "msint"){
    data_sel <- data_sel[,c("RawFile",'Proteins','ms_unique_pep', 'ms_intensity')]
    # RENAMING!
    names(data_sel) <- c('id','ms_uniprot_ac','ms_unique_pep','ms_intensity')
    # remove interactions with ms_intensity=0
    if(any(data_sel$ms_intensity == 0)) { data_sel <- data_sel[-which(data_sel$ms_intensity==0),]}
  }else{
    stop("\n << ",quant_variable," >> NOT ALLOWED. ONLY msint OR msspc ALLOWED\n")
  }

  # annotate proteins and add Masses for Mist
  species_split = unlist(strsplit(species, "-"))
  Uniprot = NULL
  for(org in species_split){
    cat(sprintf("\tLOADING %s\n",org))
    tmp <- read.delim(sprintf("%s/uniprot_protein_descriptions_%s.txt",uniprot_dir,org), stringsAsFactors= FALSE, quote="")
    if(is.null(Uniprot)){
      Uniprot = as.data.frame(tmp)
    }else{
      Uniprot = rbind(Uniprot, tmp)  
    }
  }
  cat('\tANNOTATING RESULTS...\n')
  results_annotated = merge(data_sel, Uniprot, all.x= TRUE, by.x='ms_uniprot_ac', by.y='Entry')
  
  # Adding the Mass column: it is required for MIST for preprocessing!
  # For that, we will calculate the mass of the whole protein just taking the average molecular 
  # weight of an amino acid: 110Da
  results_annotated$Mass <- results_annotated$Length*110
  
  # create output directory if it doesnt exist
  if(!dir.exists(dirname(output_file))){
    dir.create(dirname(output_file), recursive = TRUE)
  }
  
  write.table(results_annotated, file=output_file, eol='\n', sep='\t', quote= FALSE, row.names= FALSE, col.names= TRUE)
  cat('\n>> MIST FILES CREATED!\n')
  
  keysout <- subset(keys, select = c(RawFile, Condition))
  keysname <- gsub(".txt",'_mistkeys.txt', output_file)
  write.table(keysout, file=keysname, eol = '\n', sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
  cat('>> MIST keys FILE ALSO CREATED. Enjoy your day!\n\n')
}

