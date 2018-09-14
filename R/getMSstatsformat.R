# ------------------------------------------------------------------------------
#' @title Generate MSstats format object and file
#' 
#' @description Takes as input a reduced version of the Evidence file and 
#' generates the input data.frame required by MSstats. 
#' It processes fractionated data as well.
#' @param data_f (data.frame) of the filtered Maxquant evidence file.
#' @param fraction (boolean) 1 or 0 option to specified whether or not
#'  is a fractionated experiment
#' @param datafile (char) The evidence file name (to generate the output file)
#' @param funfunc (char) The function to use to aggregating the data if it is a 
#' fractionated experiment (default: `sum`)
#' @keywords internal, MSstats, format, input, fractions
#' .artms_getMSstatsFormat()
.artms_getMSstatsFormat <- function(data_f, fraction, datafile, funfunc = "sum"){
  cat("\n>> ADAPTING THE DATA TO MSSTATS FORMAT\n")
  
  data_f <- artms_changeColumnName(data_f, "Modified.sequence", "PeptideSequence")
  data_f$PeptideSequence <- gsub("_", "", data_f$PeptideSequence)
  cat("------- + Selecting Sequence Type: MaxQuant 'Modified.sequence' column\n")

  # DEAL WITH FRACTIONS FIRST (but in reality it is just checking, 
  # because it is doing a sum up of redundant features anyway)
  if( any(grepl("FractionKey", colnames(data_f))) & fraction){
    cat("------- + DEALING WITH FRACTIONS (sum up intensities per feature)\n")
    predmss <- aggregate(data = data_f, Intensity~Proteins+PeptideSequence+Charge+IsotopeLabelType+Condition+BioReplicate+Run, FUN = funfunc)
    predmss <- predmss[,c("Proteins", "PeptideSequence", "Charge", "IsotopeLabelType", "Condition", "BioReplicate", "Run", "Intensity")]
  }else{
    # If there are duplications, sum up
    predmss <- aggregate(data = data_f, Intensity~Proteins+PeptideSequence+Charge+IsotopeLabelType+Condition+BioReplicate+Run, FUN = sum)
    predmss <- predmss[,c("Proteins", "PeptideSequence", "Charge", "IsotopeLabelType", "Condition", "BioReplicate", "Run", "Intensity")]
  }
  
  # step required by MSstats to add 'NA' intensity values for those 
  # features not found in certain bioreplicates/runs
  # If this is not done, MSstats will still works, 
  # but it will generate a gigantic warning.
  # Using dcast from data.table because it has the option "sep" that allows to 
  # choose the 'collapse' character to use.
  cat("------- + Adding NA values for missing values (required by MSstats)\n")
  predmss_dc <- data.table::dcast(data = setDT(predmss), Proteins+PeptideSequence+Charge+IsotopeLabelType~Condition+BioReplicate+Run, value.var = "Intensity", fun.aggregate = sum, sep = "___")
  predmss_melt <- reshape2::melt(data = predmss_dc, id.vars = c('Proteins', 'PeptideSequence', 'Charge', 'IsotopeLabelType'), value.name = "Intensity")
  # And put back the condition, bioreplicate and run columns
  predmss_melt$Condition <- gsub("(.*)(___)(.*)(___)(.*)", "\\1", predmss_melt$variable)
  predmss_melt$BioReplicate <- gsub("(.*)(___)(.*)(___)(.*)", "\\3", predmss_melt$variable)
  predmss_melt$Run <- gsub("(.*)(___)(.*)(___)(.*)", "\\5", predmss_melt$variable)
  
  # After the data has been aggregated, then we add the columns
  predmss_melt$ProductCharge <- NA 
  predmss_melt$FragmentIon <- NA
  
  # Names required by MSstats
  predmss_melt <- artms_changeColumnName(predmss_melt, "Proteins", "ProteinName")
  predmss_melt <- artms_changeColumnName(predmss_melt, "Charge", "PrecursorCharge")
  
  # And re-sort it as msstats likes it
  dmss <- predmss_melt[,c("ProteinName", "PeptideSequence", "PrecursorCharge", "FragmentIon", "ProductCharge", "IsotopeLabelType", "Condition", "BioReplicate", "Run", "Intensity")]
  
  ## sanity check for zero's
  if(nrow(dmss[!is.na(dmss$Intensity) & dmss$Intensity == 0,]) > 0){
    dmss[!is.na(dmss$Intensity) & dmss$Intensity == 0,]$Intensity = NA
  }
  
  dmss <- as.data.frame(dmss)
  cat("------- + Write out the MSstats input file (-mss.txt)\n")
  write.table(dmss, file=gsub('.txt','-mss.txt',datafile), eol="\n", sep="\t", quote=F, row.names=F, col.names=T)
  return(dmss)  
}