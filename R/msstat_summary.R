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
#' @param evidence_file (char or data.frame) The filepath to the MaxQuant 
#' searched data (evidence) file (txt tab delimited file). 
#' Only works for the newer versions of the evidence file.
#' @param prot_group_file (char) The filepath to the MaxQuant
#' `proteinGroups.txt` file (txt tab delimited file) or data.frame
#' @param keys_file (char) The filepath to the keys file used with MSStats
#' (txt tab delimited file).
#' @param results_file (char) The filepath to the MSStats results file in t
#' he default long format (txt tab delimited file or data.frame).
#' @param return_df (data.frame) Whether or not to return the results
#' to the R environment upon completion. This is useful if this is being
#' used in an R pipeline and you want to feed the results directly into the
#' next stage of analysis via an R environment/terminal.
#' Regardless, the results will be written to file. Default = FALSE
#' @param verbose (logical) `TRUE` (default) shows function messages
#' @return (data.frame or txt file) with the summary
#' @keywords MaxQuant, evidence, MSStats, summary
#' @examples
#' # Testing warning if files are not submitted
#' test <- artmsMsstatsSummary(evidence_file = NULL,
#'                       prot_group_file = NULL,
#'                       keys_file = NULL,
#'                       results_file = NULL)
#' @export
artmsMsstatsSummary <- function(evidence_file,
                                prot_group_file,
                                keys_file,
                                results_file,
                                return_df = FALSE,
                                verbose = TRUE) {
  
  if(any(missing(evidence_file) | 
         missing(prot_group_file) |
         missing(keys_file) | 
         missing(results_file)))
    stop("Missed (one or many) required argument(s)
         Please, check the help of this function to find out more")
  
  
  # Check if passing in data or if passing in files
  if(verbose) message(">> GENERATING A GLOBAL SUMMARY")
  
  
  if(is.null(evidence_file) & 
     is.null(prot_group_file) & 
     is.null(keys_file)
     & is.null(results_file)){
    return("Files must not be NULL")
  }
  
  dat <- artmsMergeEvidenceAndKeys(evidence_file, 
                                    keys_file, 
                                    verbose = verbose)
  dat <- data.table(dat)
  pg <- .artms_checkIfFile(prot_group_file)
  pg <- data.table(pg)
  results <- .artms_checkIfFile(results_file)
  results <- data.table(results)
  
  # get SPECTRAL COUNTS
  message("--- Summarizing Spectral Counts ")
  if( any(grepl("MS.MS.count", colnames(dat))) ){
    dat <- artmsChangeColumnName(dataset = dat, "MS.MS.count", "MS.MS.Count")
  }
  # dat.sc <- data.table::dcast(data = dat,
  #                            Proteins ~ BioReplicate,
  #                            value.var = "MS.MS.Count",
  #                            max,
  #                            fill = NA_real_)
  
  dat.sc <- dat %>% 
    tidyr::pivot_wider(id_cols = Proteins, 
                       names_from = BioReplicate, 
                       values_from = MS.MS.Count, 
                       values_fn = list(MS.MS.Count = max))
  dat.sc <- as.data.frame(dat.sc)
  
    
  names(dat.sc)[-1] <- paste0(names(dat.sc)[-1], "_SC")
  
  # get INTENSITIES
  message("--- Summarizing Intensities ")
  ##LEGACY
  # dat.intensity1 <- data.table::dcast(data = dat,
  #                                    Proteins ~ BioReplicate,
  #                                    value.var = "Intensity",
  #                                    max,
  #                                    fill = NA_real_)
  
  dat.intensity <- dat %>%
    tidyr::pivot_wider(id_cols = Proteins, 
                       names_from = BioReplicate, 
                       values_from = Intensity, 
                       values_fn = list(Intensity = max))
  dat.intensity <- as.data.frame(dat.intensity)
  
  names(dat.intensity)[-1] <- paste0(names(dat.intensity)[-1], "_Intensity")
  
  # find the UNIQUE PEPTIDE columns
  message("--- Summarizing Unique Peptides ")
  idx <- grep("Peptide.counts..unique.", colnames(pg), fixed = TRUE)
  pg.uniqPep <- pg[, c(1, idx), with = FALSE]
  # # fix names to match the rest of the data and to merge smoothly
  # names(pg.uniqPep) <- gsub("Peptide counts (unique)", "", names(pg.uniqPep))
  # names(pg.uniqPep)[-1] <- paste0(names(pg.uniqPep)[-1], "_UniqPep")
  # names(pg.uniqPep)[1] <- "Proteins"
  
  names(pg.uniqPep)[grep("Peptide.counts..unique.", 
                         names(pg.uniqPep), 
                         fixed = TRUE)] <- "UniquePeptides"
  names(pg.uniqPep)[1] <- "Proteins"
  
  # convert RESULTS to WIDE format
  message(">> Converting Results to Wide format ")
  
  ##LEGACY
  # results_l1 = data.table::melt(data = results[, c("Protein", "Label", "log2FC", "adj.pvalue"), 
  #                                             with = FALSE], id.vars = c("Protein", "Label"))
  
  results_l <- results %>%
    dplyr::select(c(Protein, Label, log2FC, adj.pvalue)) %>%
    tidyr::pivot_longer(cols = -c(Protein, Label), 
                        values_to = "value", 
                        names_to = "variable")
                     
  ## then cast to get combinations of LFCV/PVAl and Label as columns
  ##LEGACY
  # results_w <- data.table::dcast(Protein ~ Label + variable,
  #                                data = results_l,
  #                                value.var = c("value"))
  
  results_w <- results_l %>%
    dplyr::mutate(Label_variable = paste(Label, variable, sep = "_")) %>%
    tidyr::pivot_wider(id_cols = Protein, 
                       names_from = Label_variable, 
                       values_from = value)
  
  
    
  names(results_w)[1] = "Proteins"
  
  # Combine them all together
  message(">> Bringing it all together ")
  results_summary = Reduce(function(...)
    merge(..., by = "Proteins", all.x = TRUE),
    list(pg.uniqPep, results_w, dat.sc, dat.intensity))
    
  names(results_summary)[grep("Proteins", names(results_summary))] = "Protein"
  
  # write out summary
  message(">> Writing out Summary. ")
  if (!is.data.frame(results_file) &
      !is.data.table(results_file)) {
    out_file <- gsub(".txt", "_summarized.txt", results_file)
    write.table(
      results_summary,
      out_file,
      quote = FALSE,
      row.names = FALSE,
      sep = "\t"
    )
  }
  
  if (return_df) {
    return(results_summary)
  }
  
  message(">> Summarization Completed ")
}
