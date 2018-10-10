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
#' @param output_file (char) The file name for the results
#' (must have the extension `.txt`). If empty, then the
#' results will be returned as an R object.
#' @param specie (char) The specie name.
#' Check `?artmsMapUniprot2Entrez`
#' @return An R object with the results and a file with the results (if the
#' output_file argument is provided). It contains averages of Intensity,
#' Retention Time, Caliberated Retention Time
#' @keywords MaxQuant, evidence, summary, intensity, retention time, caliberated
#' @examples
#' ave_int <- artms_avg_intensity_RT(evidence_file = artms_data_ph_evidence,
#'                                   specie = "human")
#' @export
artms_avg_intensity_RT <- function(evidence_file,
                                   protein_file = NULL,
                                   output_file = FALSE,
                                   specie) {
  # read in data
  cat(">> READING IN FILES...\n")
  dat <- .artms_checkIfFile(evidence_file, is.evidence = TRUE)
  suppressMessages(dat <-
                     artms_annotationUniprot(dat, "Proteins", sps = specie))
  
  if (!is.null(protein_file)) {
    cat(">> FILTERING OUT UNWANTED PROTEINS...\n")
    proteins <- .artms_checkIfFile(protein_file)
    
    # pull out only cases where the proteins appear
    query <- paste(t(proteins), collapse = "|")
    idx <- grep(query, dat$Proteins)
    dat <- dat[idx,]
  }
  
  cat(">> COMPUTING AVERAGES...\n")
  # Compute the average Intensity
  dat.avg <-
    aggregate(data = dat[, c("Protein",
                             "Gene",
                             "Modified.sequence",
                             "Charge",
                             "Intensity")], Intensity ~ .,
              mean, na.rm = TRUE)
  # Compute the average Retention Time
  dat.ret <-
    aggregate(data = dat[, c("Protein",
                             "Gene",
                             "Modified.sequence",
                             "Charge",
                             "Retention.time")], Retention.time ~
                ., mean, na.rm = TRUE)
  # Compute the average Calibrated retention time
  dat.cal <-
    aggregate(data = dat[, c("Protein",
                             "Gene",
                             "Modified.sequence",
                             "Charge",
                             "Calibrated.retention.time")],
              Calibrated.retention.time ~ .,
              mean,
              na.rm = TRUE)
  
  cat(">> MERGING RESULTS...\n")
  results <-
    merge(
      dat.avg,
      dat.ret,
      by = c("Protein", "Gene", "Modified.sequence", "Charge"),
      all.y = TRUE
    )
  results <-
    merge(
      results,
      dat.cal,
      by = c("Protein", "Gene", "Modified.sequence", "Charge"),
      all.y = TRUE
    )
  # add 'Avg' to names
  names(results)[5:7] = paste0("Avg_", names(results)[5:7])
  
  cat(">> SUMMARIZATION COMPLETE!!\n")
  
  if (output_file) {
    # write out results
    cat("--- WRITING OUT RESULTS TO ", output_file, "\n")
    write.table(
      results,
      output_file,
      quote = FALSE,
      row.names = FALSE,
      sep = "\t"
    )
    return(results)
  } else{
    return(results)
  }
}
