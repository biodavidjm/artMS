# ------------------------------------------------------------------------------
#' @title Summarize the MSStats results and data quantification
#' 
#' @description Converts the MSStats results file to wide format 
#' (unique Protein ID and columns are the comparisons), as well as adds 
#' BioReplicate information about 
#' - the Number of Unique Peptides, 
#' - Spectral Counts
#' - Intensities 
#' for each protein. 
#' In cases where there are multiple values for a Protein-BioReplicate
#' pair due to minute changes in sequence, the maximum value is taken for the 
#' pair. Any pairs without a value are assigned a value of NA.
#' @param evidence_file (char) The filepath to the MaxQuant searched data 
#' (evidence) file (txt tab delimited file). Only works in the newer versions 
#' of the evidence file
#' @param prot_group_file (char) The filepath to the MaxQuant 
#' `proteinGroups.txt`
#'  file (txt tab delimited file).
#' @param keys_file (char) The filepath to the keys file used with MSStats 
#' (txt tab delimited file).
#' @param results_file (char) The filepath to the MSStats results file in t
#' he default long format (txt tab delimited file).
#' @param return_df (data.frame) Whether or not to return the results 
#' to the R environment upon completion. This is useful if this is being 
#' used in an R pipeline and you want to feed the results directly into the 
#' next stage of analysis via an R environment/terminal. 
#' Regardless, the results will be written to file. Default = FALSE
#' @return (data.frame or txt file) with the summary
#' @keywords MaxQuant, evidence, MSStats, summary
#' @examples \donttest{
#' artms_msstats_summary(evidence_file = "FLU-THP1-H1N1-AB-evidence.txt", 
#'                 prot_group_file = "proteinGroups.txt", 
#'                 keys_file = "FLU-THP1-H1N1-AB-keys.txt", 
#'                 results_file = "results/FLU-THP1-H1N1-AB-results.txt", 
#'                 return_df = FALSE)
#' }
#' @export
artms_msstats_summary <- function(evidence_file, 
                            prot_group_file, 
                            keys_file, 
                            results_file, 
                            return_df = FALSE) {
    
    # Check if passing in data or if passing in files
    cat(">> Getting data ...\n")
    evidence <- .artms_checkIfFile(evidence_file, is.evidence = T)
    pg <- .artms_checkIfFile(prot_group_file)
    keys <- .artms_checkIfFile(keys_file)
    results <- .artms_checkIfFile(results_file)
    
    # add CONDITIONS to the evidence file
    dat <- merge(evidence, keys, by.x = "Raw file", by.y = "RawFile")
    # get SPECTRAL COUNTS
    cat(">>   Summarizing Spectral Counts\n")
    dat.sc <- data.table::dcast(data = dat, Proteins ~ BioReplicate, value.var = "MS/MS Count", max, fill = NA_real_)
    names(dat.sc)[-1] <- paste0(names(dat.sc)[-1], "_SC")
    # get INTENSITIES
    cat(">>   Summarizing Intensities\n")
    dat.intensity <- dcast(data = dat, Proteins ~ BioReplicate, value.var = "Intensity", max, fill = NA_real_)
    names(dat.intensity)[-1] <- paste0(names(dat.intensity)[-1], "_Intensity")
    
    # find the UNIQUE PEPTIDE columns
    cat(">> Summarizing Unique Peptides\n")
    idx <- grep("Peptide counts (unique)", colnames(pg), fixed = TRUE)
    pg.uniqPep <- pg[, c(1, idx), with = F]
    # # fix names to match the rest of the data and to merge smoothly
    # names(pg.uniqPep) <- gsub("Peptide counts (unique)", "", names(pg.uniqPep))
    # names(pg.uniqPep)[-1] <- paste0(names(pg.uniqPep)[-1], "_UniqPep")
    # names(pg.uniqPep)[1] <- "Proteins"
    
    names(pg.uniqPep)[grep("Peptide counts (unique)", names(pg.uniqPep), fixed = TRUE)] <- "UniquePeptides"
    names(pg.uniqPep)[1] <- "Proteins"
    
    # convert RESULTS to WIDE format
    cat(">> Converting Results to Wide format\n")
    results_l = melt(data = results[, c("Protein", "Label", "log2FC", "adj.pvalue"), with = F], id.vars = c("Protein", "Label"))
    ## then cast to get combinations of LFCV/PVAl and Label as columns
    results_w <- dcast(Protein ~ Label + variable, data = results_l, value.var = c("value"))
    names(results_w)[1] = "Proteins"
    
    # Combine them all together
    cat(">> Bringing it all together\n")
    results_summary = Reduce(function(...) merge(..., by = "Proteins", all.x = T), list(pg.uniqPep, results_w, dat.sc, dat.intensity))
    names(results_summary)[grep("Proteins", names(results_summary))] = "Protein"  # updating for annotation purposes later
    
    # write out summary
    cat(">> Writing out Summary.\n")
    if (!is.data.frame(results_file) & !is.data.table(results_file)) {
        out_file <- gsub(".txt", "_summarized.txt", results_file)
        write.table(results_summary, out_file, quote = F, row.names = F, sep = "\t")
    }
    
    if (return_df){
      return(results_summary)
    }
        
    cat(">> Summarization Complete! Have a nice day :)\n")
}




