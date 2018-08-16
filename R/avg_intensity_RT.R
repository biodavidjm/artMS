# ==============================================================================
# Input an evidence file from MaxQuant and a file containing a list of proteins 
# of interest.  
# Output a file containing averages of Intensity, Retention Time, Caliberated 
# Retention Time

#' @title Summarize average intensity and retention time per protein
#' @description This will summarize the evidence file and report back the 
#' average intensity, average retention time, and the average caliberated 
#' retention time. If a list of proteins is provided, then only those proteins 
#' will be summarized and returned.
#' @param evidence_file The filepath to the MaxQuant searched data (evidence) 
#' file (txt tab delimited file).
#' @param protein_file The filepath to a file or vector conatining a list of 
#' proteins of interest.
#' @param output_file The filepath to where you want the results saved. If no 
#' file path is given, then the results will be returned as an R object.
#' @keywords MaxQuant, evidence, summary, intensity, retention time, caliberated
#' avg_intensity_RT()
#' @export
avg_intensity_RT <- function(evidence_file, protein_file = NULL, output_file = FALSE) {
    # read in data
    cat("READING IN FILES...\n")
    dat <- checkIfFile(evidence_file, is.evidence = T)
    names(dat) <- gsub(" ", "_", names(dat))
    
    # proteins <- read.delim(protein_file, sep='\t', stringsAsFactors=F, header=F)
    if (!is.null(protein_file)) {
        cat("FILTERING OUT UNWANTED PROTEINS...\n")
        proteins <- checkIfFile(protein_file)
        
        # pull out only cases where the proteins appear
        query <- paste(t(proteins), collapse = "|")
        idx <- grep(query, dat$Proteins)
        dat <- dat[idx, ]
    }
    
    cat("COMPUTING AVERAGES...\n")
    # Compute the average Intensity
    dat.avg <- aggregate(data = dat[, c("Proteins", "Gene_names", "Modified_sequence", "Charge", "Intensity")], Intensity ~ ., 
        mean, na.rm = T)
    # Compute the average Retention Time
    dat.ret <- aggregate(data = dat[, c("Proteins", "Gene_names", "Modified_sequence", "Charge", "Retention_time")], Retention_time ~ 
        ., mean, na.rm = T)
    # Compute the average Calibrated retention time
    dat.cal <- aggregate(data = dat[, c("Proteins", "Gene_names", "Modified_sequence", "Charge", "Calibrated_retention_time")], 
        Calibrated_retention_time ~ ., mean, na.rm = T)
    
    cat("\tMERGING RESULTS...\n")
    results <- merge(dat.avg, dat.ret, by = c("Proteins", "Gene_names", "Modified_sequence", "Charge"), all.y = T)
    results <- merge(results, dat.cal, by = c("Proteins", "Gene_names", "Modified_sequence", "Charge"), all.y = T)
    # add 'Avg' to names
    names(results)[5:7] = paste0("Avg_", names(results)[5:7])
    
    cat("SUMMARIZATION COMPLETE!!\n")
    
    if (output_file) {
        # write out results
        cat("WRITING OUT RESULTS TO\n\t", gsub(".txt", "_IntRT_summary.txt", evidence_file), " ...\n")
        write.table(results, gsub(".txt", "_IntRT_summary.txt", evidence_file), quote = F, row.names = F, sep = "\t")
    } else {
        return(results)
    }
}


# # defining the structure of t spec = matrix(c( 'help' , 'h', 0, 'logical', 'available arguments (this screen)', 'file' , 'f',
# 1, 'character', 'Path to the evidence file to summarize.', 'proteins' , 'p', 1, 'character', 'Path to file containing a list
# of proteins to be searched for.'), byrow=T, ncol=5) opt = getopt(spec = spec, opt = commandArgs(TRUE), command =
# get_Rscript_filename(), usage = FALSE, debug = FALSE) # Handle the case where help is used if ( !is.null(opt$help) ) {
# cat(getopt(spec, usage=TRUE)); q(status=1); } file_summary(opt$file, opt$proteins)






