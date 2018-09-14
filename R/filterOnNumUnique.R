# ------------------------------------------------------------------------------
#' @title Filter MaxQuant data based on Number of Unique Peptides
#' 
#' @description Remove all condition-protein instance where the number of 
#' unique peptides is below the `min_pep` threshold. 
#' The results are written out to a file that is named similar to the evidence 
#' file, but also contains the filtering cutoff used.
#' @param evidence_file (char) The filepath to the MaxQuant searched data 
#' (evidence) file (txt tab delimited file).
#' @param prot_groups_file The filepath to the MaxQuant searched 
#' `proteinGroups.txt` (txt tab delimited file).
#' @param keys_file (char) The filepath to the keys file. The condition names 
#' must match the Condition names provided to the proteinGroups file replicates.
#' @param min_pep (int) The minimum number of peptides that need 
#' to be found with at least one bait-prey protein pair to consider it a 
#' real hit. Anything not meeting this criteria is removed before scoring.
#' @param verbose (Boolean) Whether or not to display status reports. 
#' Default = TRUE
#' @param return.results (Boolean) Whether or not to return the results to 
#' the R environment upon completion. This is useful if this is being used 
#' in an R pipeline and you want to feed the results directly into the next 
#' stage of analysis via an R environment/terminal. 
#' Regardless, the results will be written to file. Default = FALSE
#' @keywords MaxQuant, evidence, MSStats
#' artms_filterNumUniqPep()
#' @export
artms_filterNumUniqPep <- function(evidence_file, 
                                   keys_file, 
                                   prot_groups_file, 
                                   min_pep, 
                                   verbose = TRUE, 
                                   return.results = FALSE) {
    
    # Check if passing in data or if passing in files
    x <- .artms_checkIfFile(evidence_file, is.evidence = T)
    x.order <- names(x)
    
    pg <- .artms_checkIfFile(prot_groups_file)
    keys <- .artms_checkIfFile(keys_file)
    
    
    # Get Unique Peptides ~~~~~~~~~~~~~~~~~~~~~~
    if (verbose) 
        cat(">> Retrieving Number of Unique Peptides\n")
    # find the unique peptide columns
    idx <- grep("Unique peptides ", names(pg))
    pg.up <- pg[, c(1, idx), with = F]
    
    if (verbose) 
        cat(">> Finding condition-protein instances based on Number of Unique Peptides >= ", min_pep, "\n")
    # melt into form to filter
    pg.up <- melt(pg.up, id.vars = names(pg.up)[1], measure.vars = names(pg.up)[-1], variable.name = "Bait", value.name = "unique_peptides", 
        variable.factor = F, value.factor = F)
    # filter based on number of minimal peptides
    pg.up <- pg.up[pg.up$unique_peptides >= min_pep, ]
    # fix the bait names to *hopefully* match the keys file. This may change as versions of MQ change.
    pg.up <- pg.up[, `:=`(Bait, gsub("Unique peptides ", "", Bait))]
    
    # keep a list of all the bait-prey pairs that pass
    pg.pairs <- unique(pg.up[, 1:2])
    pg.pairs <- pg.pairs[, `:=`(pairs, paste(`Protein IDs`, Bait, sep = "-"))]
    
    # Evidence and Keys ~~~~~~~~~~~~~~~~~~~~
    if (verbose) 
        cat(">> Filtering based on Number of Unique Peptides\n")
    # Fix the condition labels. There is no 'RNase' distinction in the pg file. There may also be others as different datasets
    # come.
    keys$BioReplicate <- gsub("_RNase", "", keys$BioReplicate)
    # merge evidence and keys files for proper bait-prey pair naming
    x <- merge(x, keys[, c("Raw.file", "BioReplicate")], by.x = "Raw file", by.y = "Raw.file")
    # create a paired name to compare with
    x$pairs <- paste(x$Proteins, x$BioReplicate, sep = "-")
    
    
    # keep only the data where the pairs make the threshold
    results <- x[x$pairs %in% pg.pairs$pairs, ]
    
    # Preserver the order and column names of the original evidence file
    results <- results[, x.order, with = F]
    
    # write out the results
    if (verbose) 
        cat(">> Writing out results\n")
    out_file <- gsub(".txt", paste0("_filtered_min", min_pep, "uniqPep.txt"), evidence_file)
    write.table(results, out_file, quote = F, row.names = F, sep = "\t")
    
    if (return.results) 
        return(results)
}



