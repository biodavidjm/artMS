# ------------------------------------------------------------------------------
#' @title MaxQuant evidence file to SAINTexpress format
#' 
#' @description Converts the MaxQuant evidence file to the 3 required files 
#' by SAINTexpress. One can choose to either use the `spectral counts`
#' (use `msspc`) or the `intensities` (use `msint`) for the analysis. 
#' @param input_file (char) Maxquant evidence file
#' @param keys_file (char) Keys file with a SAINT column specifying 
#' test (`T`) and control (`C`) conditions
#' @param output_file (char) Output file name
#' @param ref_proteome_file (char) Reference proteome path file name in 
#' fasta format
#' @param quant_variable (char) choose either 
#' - `msspc` (spectral counts, default) or
#' - `msint` (MS Intensity)
#' @return The 3 required files by SAINTexpress
#' @keywords SAINT, SAINTexpress, APMS
#' @examples \donttest{
#' artms_evidenceToSaintExpressFormat(input_file = "a549-PB1-evidence.txt", 
#'              keys_file = "a549-PB1-keys.txt", 
#'              ref_proteome_file = "fluomics-uniprot-hsa_20170516.fasta", 
#'              quant_variable = "msint", 
#'              output_file = "a549-PB1-saintexpress.txt")
#' }
#' @export
artms_evidenceToSaintExpressFormat <- function(input_file, 
                                               keys_file, 
                                               ref_proteome_file, 
                                               quant_variable='msspc', 
                                               output_file){
  
  cat(">> CONVERTING TO SAINTexpress FORMAT\n")
  
  data = fread(input_file, integer64 = 'double')
  keys = fread(keys_file, integer64 = 'double')
  
  ## write baits in format
  ## hIP101-10       PV_2C_co_uni    T
  
  saint_baits = keys[,c('BioReplicate','Condition','SAINT'),with= FALSE]
  
  ## write interactions in format
  ## hIP101-10       PV_2C_co_uni    Q9NTG7  1
  
  tryCatch(setnames(data, 'Raw file', 'RawFile'), error=function(e) cat("--- 'Raw file' not found\n"))
  tryCatch(setnames(keys, 'Raw.file', 'RawFile'), error=function(e) cat("--- 'Raw.file' not found\n"))
  
  cat('>> VERIFYING DATA AND KEYS\n')
  if(any(!c('RawFile','IsotopeLabelType','Condition','BioReplicate','Run','SAINT') %in% colnames(keys))){
    stop('COLNAMES IN KEYS NOT CONFORM TO SCHEMA\n\tRawFile\tIsotopeLabelType\tCondition\tBioReplicate\tRun\tSAINT\n')
  } 
  if(!'IsotopeLabelType' %in% colnames(data)) data[,IsotopeLabelType:='L']
  data <- artms_mergeEvidenceAndKeys(data, keys, by = c('RawFile','IsotopeLabelType'))
  data_f <- artms_filterEvidenceContaminants(data)
  data_f <- .artms_removeMaxQProteinGroups(data_f) ## do we want this or not?
  
  cat(">> AGGREGATING ON", quant_variable,"VALUES...\n")
  ## aggregate over technical replicates if necessary
  if(quant_variable == 'msspc'){
    setnames(data_f, 'MS/MS Count', 'spectral_counts')
    data_f_agg <- aggregate(spectral_counts ~ BioReplicate+Condition+Proteins+Sequence+Charge,data=data_f, FUN = max)
    data_f_agg <- aggregate(spectral_counts ~ BioReplicate+Condition+Proteins,data=data_f_agg, FUN = sum)
  }else if(quant_variable=='msint'){
    data_f_agg <- aggregate(Intensity ~ BioReplicate+Condition+Proteins+Sequence+Charge, data=data_f, FUN = max)
    data_f_agg <- aggregate(Intensity ~ BioReplicate+Condition+Proteins, data=data_f_agg,FUN = sum)
  }else{
    stop("\nERROR!! Wrong value for variable to quantify. Please use 'msspc' or 'msint'")
  }
  
  ## IP name, bait name, prey name, and spectral counts or intensity values
  saint_interactions = data_f_agg
  
  ## write preys in format
  ## Q9NTG7  43573.5 Q9NTG7
  
  ref_proteome = read.fasta(file = ref_proteome_file, 
                            seqtype = "AA", 
                            as.string = TRUE,
                            set.attributes = TRUE, 
                            legacy.mode = TRUE, 
                            seqonly = FALSE, 
                            strip.desc = FALSE)
  p_lengths = c()
  p_names = c()
  for(e in ref_proteome){
    p_lengths = c(p_lengths, nchar(e[1]))
    p_names = c(p_names, attr(e,'name'))
  }
  ref_table = data.table(names=p_names, lengths=p_lengths)
  ref_table[,uniprot_ac:=gsub('([a-z,0-9,A-Z]+\\|{1})([A-Z,0-9,\\_]+)(\\|[A-Z,a-z,0-9,_]+)','\\2',names)]
  ref_table[,uniprot_id:=gsub('([a-z,0-9,A-Z]+\\|{1})([a-z,0-9,A-Z]+\\|{1})([A-Z,a-z,0-9,_]+)','\\3',names)]
  
  unique_preys = data.table(uniprot_ac=unique(data_f_agg$Proteins))
  saint_preys = ref_table[,c('uniprot_ac','lengths','uniprot_id'),with= FALSE]
  saint_preys = merge(unique_preys, saint_preys, by='uniprot_ac', all.x= TRUE)
  missing_lengths = nrow(saint_preys[is.na(saint_preys$uniprot_id),])
  saint_preys[is.na(saint_preys$uniprot_id),]$uniprot_id=saint_preys[is.na(saint_preys$uniprot_id),]$uniprot_ac  
  if(missing_lengths>0){
    cat(sprintf("--- WARNING! COMPUTING %s MISSING LENGTHS WITH THE MEDIAN LENGTH FROM THE DATASET\n",missing_lengths))
    saint_preys[is.na(saint_preys$lengths),]$lengths=median(saint_preys$lengths,na.rm = TRUE)
  }
  
  # Check if output directory exists and create it if not.
  if(!dir.exists(dirname(output_file))){
    dir.create(dirname(output_file), recursive= TRUE)
  }
  
  ## WRITE
  write.table(saint_baits, file = gsub('.txt','-saint-baits.txt',output_file), eol='\n',sep='\t', quote= FALSE, row.names= FALSE, col.names= FALSE)
  write.table(saint_preys, file = gsub('.txt','-saint-preys.txt',output_file), eol='\n',sep='\t', quote= FALSE, row.names= FALSE, col.names= FALSE)
  write.table(saint_interactions, file = gsub('.txt','-saint-interactions.txt',output_file), eol='\n',sep='\t', quote= FALSE, row.names= FALSE, col.names= FALSE)
  cat(">> OUTPUT FILES:\n")
  cat("--- ", gsub('.txt','-saint-baits.txt', output_file), "\n")
  cat("--- ", gsub('.txt','-saint-preys.txt', output_file), "\n")
  cat("--- ", gsub('.txt','-saint-interactions.txt', output_file, "\n"))
}


