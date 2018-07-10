#' @import data.table
#' @title Check if an input is a file or a data object
#' @description This function is used in order to make it so a user can submit
#'  either a path to a data file or a data object in data.frame or data.table 
#'  form.
#' @param input_file The filepath/object to be checked.
#' @param is.evidence Whether or not the file to be read in is an 
#' evidence file. This will assign proper classes to the evidence file 
#' when being read in.
#' @keywords file, evidence, input
#' checkIfFile()
#' @export
checkIfFile <- function(input_file, is.evidence = FALSE) {
    # check if already a data.frame or data.table
    if (is.data.table(input_file) | is.data.frame(input_file)) {
        x <- data.table(input_file)
    } else if (tryCatch(file.exists(input_file), 
                        error = function(e) return(FALSE))) {
        # check if file path is legit
        if (is.evidence) {
            x <- read_evidence_file(input_file)
        } else {
            x <- fread(input_file, integer64 = "double")
        }
    } else {
        stop("There's something wrong with the file/object you submitted:\n\t", 
             input_file, "\nPlease check that the file directory is correct, 
             or that the data is in `data.frame` or `data.table` format.")
    }
    return(x)
}

#' @title Check MaxQuant's evidence file version
#' @description MaxQuant introduced changes in the column names and number
#' of columns for the evidence file. This function check for the MaxQuant 
#' version
#' @param evidence.txt the evidence file name
#' @keywords file, evidence, input, check version
#' checkEvidenceVersion()
#' @export
checkEvidenceVersion <- function(evidence_file){
  evidence <-read.delim(evidence_file, sep = "\t", quote = "", header = T, stringsAsFactors = F)
  #old file formats
  oldfile <-
    c(
      "Sequence",
      "Length",
      "Modifications",
      "Modified.sequence",
      "Oxidation..M..Probabilities",
      "Oxidation..M..Score.Diffs",
      "Acetyl..Protein.N.term.",
      "Oxidation..M.",
      "Proteins",
      "Leading.Proteins",
      "Leading.Razor.Protein",
      "Gene.Names",
      "Protein.Names",
      "Fasta.headers",
      "Type", 
      "Raw.file",
      "MS.MS.m.z",
      "Charge",
      "m.z",
      "Mass",
      "Resolution", 
      "Uncalibrated...Calibrated.m.z..ppm.",
      "Mass.Error..ppm.",
      "Uncalibrated.Mass.Error..ppm.",
      "Max.intensity.m.z.0",
      "Retention.time",
      "B.mixture",
      "Retention.Length",
      "Calibrated.Retention.Time",
      "Calibrated.Retention.Time.Start",
      "Calibrated.Retention.Time.Finish",
      "Retention.Time.Calibration",
      "Match.Time.Difference",
      "PIF",
      "Fraction.of.total.spectrum",
      "Base.peak.fraction",
      "PEP",
      "MS.MS.Count",
      "MS.MS.Scan.Number", 
      "Score",
      "Delta.score",
      "Combinatorics",
      "Intensity",
      "Reverse",
      "Contaminant",
      "id",
      "Protein.group.IDs",
      "Peptide.ID",
      "Mod..peptide.ID",
      "MS.MS.IDs",
      "Best.MS.MS",
      "AIF.MS.MS.IDs",
      "Oxidation..M..site.IDs"
    )
  #new file format
  newfile <- c(
    "Sequence",
    "Length",
    "Modifications",
    "Modified.sequence",
    "Oxidation..M..Probabilities",
    "Oxidation..M..Score.Diffs",
    "Acetyl..Protein.N.term.",
    "Oxidation..M.",
    "Missed.cleavages",
    "Proteins",
    "Leading.proteins",
    "Leading.razor.protein", 
    "Gene.names",
    "Protein.names",
    "Type",
    "Raw.file",
    "MS.MS.m.z",
    "Charge",
    "m.z",
    "Mass",
    "Resolution", 
    "Uncalibrated...Calibrated.m.z..ppm.",
    "Uncalibrated...Calibrated.m.z..Da.",
    "Mass.Error..ppm.",
    "Mass.Error..Da.",
    "Uncalibrated.Mass.Error..ppm.",
    "Uncalibrated.Mass.Error..Da." ,
    "Max.intensity.m.z.0",
    "Retention.time",
    "Retention.length",
    "Calibrated.retention.time",
    "Calibrated.retention.time.start",
    "Calibrated.retention.time.finish",
    "Retention.time.calibration",
    "Match.time.difference",
    "Match.m.z.difference",
    "Match.q.value",
    "Match.score",
    "Number.of.data.points",
    "Number.of.scans",
    "Number.of.isotopic.peaks",
    "PIF",
    "Fraction.of.total.spectrum",
    "Base.peak.fraction",
    "PEP",
    "MS.MS.Count",
    "MS.MS.Scan.Number",
    "Score",
    "Delta.score",
    "Combinatorics",
    "Intensity",
    "Reverse",
    "Potential.contaminant",
    "id",
    "Protein.group.IDs",
    "Peptide.ID",
    "Mod..peptide.ID",
    "MS.MS.IDs",
    "Best.MS.MS",
    "AIF.MS.MS.IDs",
    "Oxidation..M..site.IDs"
  )
  cat("+--- File opened successfully\n")
  #testing file type
  if (all(c(colnames(rawEvidence)==c(newfile)))==T){
    cat("+--- Newer version of MaxQuant evidence file\n")
  }else if(all(c(colnames(rawEvidence)==c(oldfile)))==T){
    cat("+--- Old version of MaxQuant evidence file\n")
  }else{
    stop("(-) File type unknown!\n")
  }
}