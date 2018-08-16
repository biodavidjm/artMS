#' @title MaxQuant PTM to MSStats.
#' @description Converts a MaxQuant PTM Sites file into a MSStats compatible file.
#' @param ptmsites_file The filepath/object to the PTM Sites.
#' @param is.evidence Whether or not the file to be read in is an evidence file. This will assign proper classes to the evidence file when being read in.
#' @keywords file, evidence, input, sites, PTM, post translational modification, STY
#' mq2mss_PTM()
#' @export
mq2mss_PTM = function(ptmsites_file) {
    
    cat("READING IN DATA...")
    if (file.exists(ptmsites_file)) {
        x = read.delim(ptmsites_file, stringsAsFactors = F)
    } else if (is.data.frame(ptmsites_file)) {
        x = ptmsites_file
    } else {
        stop("There was a problem accessing the ptmsites_file. Please make sure it is either a correct file path or a data.frame object.")
    }
    cat("DONE\n")
    
    
    cat("WRANGLING DATA...\n")
    cat("\tCREATING SITE INFORMATION...\n")
    # explode the Intensities to different lines
    x.long = melt(x[, c(grep("^Phospho\\.\\.STY\\.\\.Probabilities$|Intensity\\.[A-Za-z0-9]+__|^id$|^Charge$", names(x)))], id = c("Phospho..STY..Probabilities", 
        "id", "Charge"), value.name = "Intensity", variable.name = "Condition")
    # add in site locations to the identifier so we can map everything back later
    x.long$ProteinName = paste0(x$Phospho..STY..Probabilities, "_", x$id)
    
    cat("\tORGANIZING COLUMNS...\n")
    # remove 0 values
    x.long = x.long[x.long$Intensity > 0, ]
    x.long = unique(x.long)
    
    # simplify the condition name
    x.long$Condition = gsub("Intensity\\.", "", x.long$Condition)
    
    # Get BioReplicate
    x.long$BioReplicate = gsub("(.*)(___.*)", "\\1", x.long$Condition)  #   <-------- CHANGE WHEN WE HAVE AN ESTABLISHED STYLE
    # Create run number - The reun information was already aggregated in the SITEs file, so we don't need to worry about that here
    x.long$Run = 1  #x.long$id  --  too many id's
    # Clean up condition
    x.long$Condition = gsub("([A-Za-z0-9]+)([0-9]+)(___)(.*)", "\\1", x.long$Condition)  #   <-------- CHANGE WHEN WE HAVE AN ESTABLISHED STYLE
    
    x.long$PeptideSequence = x.long$Phospho..STY..Probabilities
    x.long$PrecursorCharge = x.long$FragmentIon = NA
    x.long$ProductCharge = x.long$Charge
    x.long$IsotopeLabelType = "L"
    
    cat("\tTRIMMING DATA...\n")
    # Save only the necessary columns
    mss = x.long[, c("ProteinName", "PeptideSequence", "PrecursorCharge", "FragmentIon", "ProductCharge", "IsotopeLabelType", "Condition", 
        "BioReplicate", "Run", "Intensity")]
    
    cat("WRITING OUT MSS FILE...\n")
    # Write out the results
    write.table(mss, gsub(".txt", "_mss.txt", ptmsites_file), quote = F, row.names = F)
    
    cat("MSS COMPATIBLE FILE SAVED. \n\tHAVE A NICE DAY :)")
    return(mss)
}



