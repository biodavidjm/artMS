# ------------------------------------------------------------------------------
# @title Generate MSstats format object and file
#
# @description Takes as input a reduced version of the Evidence file and
# generates the input data.frame required by MSstats.
# It processes fractionated data as well.
# @param data_f (data.frame) of the filtered Maxquant evidence file.
# @param output_name (char) Output file name (to generate the output files). 
# '.txt' extension required
# @param data_object (logical) if TRUE the output_name cannot be the evidence file
# since it would be an data_object
# @param verbose (logical) `TRUE` (default) shows function messages
# @return (data.frame) MSstats compatible format
# @keywords internal, MSstats, format, input, fractions
.artms_getMSstatsFormat <- function(data_f, 
                                    output_name, 
                                    data_object = FALSE,
                                    verbose = TRUE) {
  
  Run = PeptideSequence = Condition_BioReplicate_Run = NULL
  
  if(verbose) message(">> CONVERTING THE DATA TO MSSTATS FORMAT ")

  if(any(missing(data_f) | 
         missing(output_name)))
    stop("Missed (one or many) required argument(s)
         Please, check the help of this function to find out more")
    
  if(data_object){
    output_name <- "artms_evidence.txt"
  }else{
    if (!grepl(".txt", output_name)) {
      stop("Argument <output_file> must have the extension '.txt'")
    }
  }
  
  data_f <- artmsChangeColumnName(data_f, 
                                  "Modified.sequence", 
                                  "PeptideSequence")
  
  data_f$PeptideSequence <- gsub("_", "", data_f$PeptideSequence)
  
  if(verbose) 
    message("-- Selecting Sequence Type: MaxQuant 'Modified.sequence' column")
  
  # DEAL WITH FRACTIONS FIRST 
  if( !("Fraction" %in% colnames(data_f)) ){
    if(verbose) message("\t(+) <Fraction> column added (with value 1, MSstats requirement)")
    data_f$Fraction <- 1
  }
  
  predmss <- data_f[, c("Proteins",
                        "PeptideSequence",
                        "Charge",
                        "IsotopeLabelType",
                        "Condition",
                        "BioReplicate",
                        "Fraction",
                        "Run",
                        "Intensity")]
  
  # step required by MSstats to add 'NA' intensity values for those
  # features not found in certain bioreplicates/runs
  # If this is not done, MSstats will still works,
  # but it will generate a gigantic warning.
  if(verbose)
    message("-- Adding NA values for missing values (required by MSstats) ")
  
  ##LEGACY
  # predmss_dc <- data.table::dcast(data = setDT(predmss),
  #                                 Proteins+PeptideSequence+Charge+IsotopeLabelType~Condition+
  #                                   BioReplicate + Run,
  #                                 value.var = "Intensity",
  #                                 fun.aggregate = sum,
  #                                 sep = "___")
  predmss_dc <- predmss %>% 
    dplyr::mutate(Condition_BioReplicate_Run = paste(Condition, BioReplicate, Run, Fraction, sep = "___") ) %>%
    tidyr::pivot_wider(id_cols = c(Proteins, PeptideSequence, Charge, IsotopeLabelType), 
                       names_from = Condition_BioReplicate_Run, 
                       values_from = Intensity, 
                       values_fn = list(Intensity = sum), values_fill = list(Intensity = NA))
                    
  
  ##LEGACY
  # predmss_melt <- data.table::melt(data = predmss_dc,
  #                                  id.vars = c('Proteins', 
  #                                              'PeptideSequence', 
  #                                              'Charge', 
  #                                              'IsotopeLabelType'),
  #                                  value.name = "Intensity")
  predmss_melt <- predmss_dc %>%
    tidyr::pivot_longer(cols = -c(Proteins, PeptideSequence, Charge, IsotopeLabelType), 
                        names_to = "variable", 
                        values_to = "Intensity")
  
  # And put back the condition, bioreplicate and run columns
  predmss_melt$Condition <- gsub("(.*)(___)(.*)(___)(.*)(___)(.*)", "\\1", predmss_melt$variable)
  predmss_melt$BioReplicate <- gsub("(.*)(___)(.*)(___)(.*)(___)(.*)", "\\3", predmss_melt$variable)
  predmss_melt$Run <- gsub("(.*)(___)(.*)(___)(.*)(___)(.*)", "\\5", predmss_melt$variable)
  predmss_melt$Fraction <- gsub("(.*)(___)(.*)(___)(.*)(___)(.*)", "\\7", predmss_melt$variable)
  
  # After the data has been aggregated, then we add the columns
  predmss_melt$ProductCharge <- NA
  predmss_melt$FragmentIon <- NA
  
  # Names required by MSstats
  predmss_melt <- artmsChangeColumnName(predmss_melt, "Proteins", "ProteinName")
  predmss_melt <- artmsChangeColumnName(predmss_melt, "Charge", "PrecursorCharge")

  # And re-sort it as msstats likes it
  dmss <- predmss_melt[, c("ProteinName",
                           "PeptideSequence",
                           "PrecursorCharge",
                           "FragmentIon",
                           "ProductCharge",
                           "IsotopeLabelType",
                           "Condition",
                           "BioReplicate",
                           "Run",
                           "Fraction",
                           "Intensity")]
  
  ## sanity check for zero's
  if ( nrow(dmss[!is.na(dmss$Intensity) & dmss$Intensity == 0, ]) > 0) {
    dmss[!is.na(dmss$Intensity) & dmss$Intensity == 0, ]$Intensity = NA
  }
  
  dmss <- as.data.frame(dmss)
  if(verbose) message("-- Write out the MSstats input file (-mss.txt) ")
  write.table(dmss,
              file = gsub('.txt', '-mss.txt', output_name),
              eol = "\n",
              sep = "\t",
              quote = FALSE,
              row.names = FALSE,
              col.names = TRUE)
    
  return(dmss)
}