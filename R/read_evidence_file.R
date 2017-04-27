
#' @title Read in Evidence File
#' @description Read in a MaxQuant searched Evidence file using data.table. This function propperly classes each column and so fread doesn't have to guess.
#' @param evidence_file The filepath to the MaxQuant searched data (evidence) file (txt tab delimited file).
#' @keywords MaxQuant, evidence
#' read_evidence_file()
#' @export
read_evidence_file <- function(evidence_file){
  cat("Reading in evidence file...\n")
  # read in the first line to get the header names
  cols <- readLines(evidence_file, 1)
  cols <- data.frame( V1 = unlist(strsplit(cols, '\t')), stringsAsFactors = F)
  cols$idx <- 1:dim(cols)[1]

  # get data frame of pre-recorded column names and their respective classes
  col.classes <- as.data.frame( matrix(c("Sequence","character","Length","integer","Modifications","character","Modified sequence","character","Oxidation (M) Probabilities","character","Oxidation (M) Score Diffs","character","Acetyl (Protein N-term)","integer","Oxidation (M)","integer","Missed cleavages","integer","Proteins","character","Leading proteins","character","Leading razor protein","character","Gene names","character","Protein names","character","Type","character","Raw file","character","Experiment","character","MS/MS m/z","numeric","Charge","integer","m/z","numeric","Mass","numeric","Resolution","numeric","Uncalibrated - Calibrated m/z [ppm]","numeric","Uncalibrated - Calibrated m/z [Da]","numeric","Mass Error [ppm]","numeric","Mass error [ppm]","numeric","Mass Error [Da]","numeric","Mass error [Da]","numeric","Uncalibrated Mass Error [ppm]","numeric","Uncalibrated mass error [ppm]","numeric","Uncalibrated Mass Error [Da]","numeric","Uncalibrated mass error [Da]","numeric","Max intensity m/z 0","numeric","Retention time","numeric","Retention length","numeric","Calibrated retention time","numeric","Calibrated retention time start","numeric","Calibrated retention time finish","numeric","Retention time calibration","numeric","Match time difference","numeric","Match m/z difference","numeric","Match q-value","numeric","Match score","numeric","Number of data points","integer","Number of scans","integer","Number of isotopic peaks","integer","PIF","numeric","Fraction of total spectrum","numeric","Base peak fraction","numeric","PEP","numeric","MS/MS Count","integer","MS/MS count","integer","MS/MS Scan Number","integer","MS/MS scan number","integer","Score","numeric","Delta score","numeric","Combinatorics","integer","Intensity","numeric","Reverse","character","Potential contaminant","character","id","integer","Protein group IDs","character","Peptide ID","integer","Mod. peptide ID","integer","MS/MS IDs","character","Best MS/MS","integer","AIF MS/MS IDs","logical","Oxidation (M) site IDs","character"          ,"Acetyl (K) Probabilities","character","GlyGly (K) Probabilities","character","Phospho (STY) Probabilities","Character","Acetyl (K) Score Diffs","character","GlyGly (K) Score Diffs","character","Phospho (STY) Score Diffs","character","Acetyl (K)","integer","GlyGly (K)","integer","Phospho (STY)","integer","Acetyl (K) site IDs","character","GlyGly (K) site IDs","character","Phospho (STY) site IDs","character"), ncol=2, byrow=T), stringsAsFactors = F)
  # merge the classes to the columns
  cols.matched = merge(cols, col.classes, by="V1", all.x=T)
  # re-order things to match the initial order
  cols.matched <- cols.matched[order(cols.matched$idx),]

  # Stop if there is an issue
  if(length(which(is.na(cols.matched$V2)))>0){
    stop(paste0("OH NO!! YOUR EVIDENCE FILE CONTAINS A COLUMN THAT I DON'T RECOGNIZE :( PLEASE TELL THE 'col.classes' IN THE read_evidence_file' FUNCTION AND ADD IN THIS NEW COLUMN(S) CALLED \n\t", paste(cols.matched$V1[which(is.na(cols.matched$V2))], collapse="\n\t"), "\n" ) )
  }

  # read in the evidence file with their classes
  x <- fread(evidence_file, integer64 = 'double', colClasses = cols.matched$V2)
  return(x)
}

