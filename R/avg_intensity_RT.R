# ==============================================================================
#' @title Summarize average intensity and retention time per protein
#' 
#' @description Input an evidence file from MaxQuant and a file 
#' containing a list of proteins of interest (optional). 
#' The function will summarize from the evidence file and report back the 
#' average intensity, average retention time, and the average caliberated 
#' retention time. If a list of proteins is provided, then only those proteins 
#' will be summarized and returned.
#' @param evidence_file (char) The filepath to the MaxQuant searched data 
#' (evidence) file (txt tab delimited file).
#' @param protein_file (char) The filepath to a file or vector conatining 
#' a list of proteins of interest.
#' @param output_file (char) The filepath to where you want the results saved 
#' (must have the extension `.txt`). If no file path is given, then the 
#' results will be returned as an R object.
#' @return An R object with the results and a file with the results (if the
#' output_file argument is provided). It contains averages of Intensity, 
#' Retention Time, Caliberated Retention Time
#' @keywords MaxQuant, evidence, summary, intensity, retention time, caliberated
#' @examples \donttest{
#' artms_avg_intensity_RT(evidence_file = "/path/to/the/evidence.txt")
#' }
#' @export
artms_avg_intensity_RT <- function(evidence_file, protein_file = NULL, output_file = FALSE) {
    # read in data
    cat(">> READING IN FILES...\n")
    dat <- .artms_checkIfFile(evidence_file, is.evidence = T)
    names(dat) <- gsub(" ", "_", names(dat))
    
    # proteins <- read.delim(protein_file, sep='\t', stringsAsFactors=F, header=F)
    if (!is.null(protein_file)) {
        cat(">> FILTERING OUT UNWANTED PROTEINS...\n")
        proteins <- .artms_checkIfFile(protein_file)
        
        # pull out only cases where the proteins appear
        query <- paste(t(proteins), collapse = "|")
        idx <- grep(query, dat$Proteins)
        dat <- dat[idx, ]
    }
    
    cat(">> COMPUTING AVERAGES...\n")
    # Compute the average Intensity
    dat.avg <- aggregate(data = dat[, c("Proteins", "Gene_names", "Modified_sequence", "Charge", "Intensity")], Intensity ~ ., 
        mean, na.rm = T)
    # Compute the average Retention Time
    dat.ret <- aggregate(data = dat[, c("Proteins", "Gene_names", "Modified_sequence", "Charge", "Retention_time")], Retention_time ~ 
        ., mean, na.rm = T)
    # Compute the average Calibrated retention time
    dat.cal <- aggregate(data = dat[, c("Proteins", "Gene_names", "Modified_sequence", "Charge", "Calibrated_retention_time")], 
        Calibrated_retention_time ~ ., mean, na.rm = T)
    
    cat(">> MERGING RESULTS...\n")
    results <- merge(dat.avg, dat.ret, by = c("Proteins", "Gene_names", "Modified_sequence", "Charge"), all.y = T)
    results <- merge(results, dat.cal, by = c("Proteins", "Gene_names", "Modified_sequence", "Charge"), all.y = T)
    # add 'Avg' to names
    names(results)[5:7] = paste0("Avg_", names(results)[5:7])
    
    cat(">> SUMMARIZATION COMPLETE!!\n")
    
    if (output_file) {
        # write out results
        outputFileNameFinal <- gsub(".txt", "_IntRT_summary.txt", evidence_file)
        cat("--- WRITING OUT RESULTS TO ", outputFileNameFinal, "\n")
        write.table(results, outputFileNameFinal, quote = F, row.names = F, sep = "\t")
        return(results)
    } else {
        return(results)
    }
}




