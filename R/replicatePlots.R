# ------------------------------------------------------------------------------
#' @title Create the replicate plots based on the pairings from
#' the replicate plot file
#'
#' @description Outputs a replicate plots based on a user provied file
#' containing the replicates to be compared. Values are based on the log2
#' value of the maximum intensities per modified sequence.
#' The "replicate file" should describe which replicates of which conditions
#' should be compared against each other. Each row represents a replicate
#' plot to be created. The file should be structured using the following
#' format and column names:
#'
#' ```
#' **condition1**|**rep1\_1**|**rep1\_2**|**condition2**|**rep2\_1**|**rep2\_2**
#'-----|-----|-----|-----|-----|-----
#' Cal33|Cal33-1|Cal33-2|HSC6|HSC6-1|HSC6-2
#' Cal33|Cal33-3|Cal33-4|HSC6|HSC6-3|HSC6-4
#' ```
#'
#' @param input_file (char) MaxQuant evidence file and location
#' @param keys_file (char) Keys file with the experimental details
#' @param replicate_file (char) Replicate file. Check Vignette for examples
#' @param out_file (char) Output .TEXT file with the intensity values of every
#' feature
#' @param prot_exp (char) Proteomics experiment. 4 options available:
#' - `AB`: (default) protein abundance
#' - `APMS`: affinity purification mass spectrometry
#' - `PH`: protein phosphorylation
#' - `UB`: protein ubiquitination (aka ubiquitylation)
#' @param verbose (logical) `TRUE` (default) shows function messages
#' @return The output file of the summary of features and intensity values
#' @keywords evidence, replica, plots
#' @examples
#' # First, let's make the "replicate file" (in a data.frame)
#' x_names <- c("condition1", "rep1_1", "rep1_2", "condition2", "rep2_1", 
#' "rep2_2")
#' x_values <- c("Cal33", "Cal33-1", "Cal33-4", "HSC6", "HSC6-2", "HSC6-3")
#' replica_info <- data.frame(t(x_values))
#' colnames(replica_info) <- x_names
#'
#' # Now let's make the plots.
#' artms_replicatePlots(input_file = artms_data_ph_evidence,
#'                      keys_file = artms_data_ph_keys,
#'                      replicate_file = replica_info,
#'                      out_file = NULL,
#'                      prot_exp = "PH")
#'
#' # Remember that if you want to see the txt results and pdf file, just
#' # change out_file = NULL" to out_file = 'output_file.pdf'
#' @export
artms_replicatePlots <- function(input_file,
                                 keys_file,
                                 replicate_file,
                                 out_file,
                                 prot_exp  = c("AB", "PH", "UB", "APMS"),
                                 verbose = TRUE) {
  if(verbose) cat(">> GENERATING CUSTOMIZED REPLICATE PLOTS\n")
  
  if(any(missing(input_file) | 
         missing(keys_file) |
         missing(replicate_file) | 
         missing(out_file)))
    stop("Missed (one or many) required argument(s)
         Please, check the help of this function to find out more")
  
  # FILTER BY PROTEOMICS EXPERIMENT
  prot_exp <- toupper(prot_exp)
  prot_exp <- match.arg(prot_exp)
  supportedExperiments <- c('AB', 'PH', 'UB', 'APMS')
  
  if (any(!prot_exp %in% supportedExperiments)) {
    stop(prot_exp, " is currently not supported.
         The experiments supported are:\n",
         sprintf('\t%s\n', supportedExperiments))
  }
  
  if (!is.null(out_file)) {
    if (!grepl(".txt", out_file)) {
      stop("<out_file> must have extension '.txt'
        Change out_file extension and try again\n"
      )
    }
  }

  if(verbose) cat("--- READING IN FILES...\n")
  # read in data
  dat <- .artms_checkIfFile(input_file)
  # keys
  keys <- .artms_checkIfFile(keys_file)
  
  # profile plot list
  repplot <- .artms_checkIfFile(replicate_file)
  
  # remove negatives from MaxQuant
  if (length(grep("__", dat$Proteins)) > 0)
    dat <- dat[-grep("__", dat$Proteins), ]
  
  # remove blank protein names
  if (any(dat$Proteins == "")) {
    dat <- dat[-which(dat$Proteins == ""), ]
  }
  
  if (prot_exp == "UB") {
    dat <- dat[grep("(gl)", dat$Modified.sequence), ]
    if(verbose) cat("--- Selecting only UB modified peptides\n")
  } else if (prot_exp == "PH") {
    dat <- dat[grep("(ph)", dat$Modified.sequence), ]
    if(verbose) cat("--- Selecting only PH modified peptides\n")
  } else if (prot_exp == "AC") {
    dat <- dat[grep("K\\(ac\\)", dat$Modified.sequence), ]
    if(verbose) cat("--- Selecting only AC modified peptides\n")
  } else if (prot_exp == "AB" | prot_exp == "APMS") {
    if(verbose) cat("--- No filtering of modified peptides\n")
  } else{
    stop(
      "\n!!! THE prot_exp IS NOT RECOGNIZED. CHECK ?artms_replicatePlots 
      TO FIND OUT THE AVAILABLE OPTIONS\n"
    )
  }
  
  # NOTE: dimensions between x and dat may differ if there is data in dat that
  # isn't in the keys file
  names(dat)[grep("Raw.file", names(dat))] <- 'RawFile'
  x <-
    merge(dat, keys[, c('RawFile', 'Condition', 'BioReplicate')], 
          by = c('RawFile'))
  
  # Put into a data matrix format
  x <-
    data.table::dcast(
      data = x,
      Proteins + Modified.sequence + Charge ~ Condition + BioReplicate,
      value.var = "Intensity",
      max,
      na.rm = TRUE
    )
  # remove cases where -Inf  is introduced
  x[x == -Inf] <- 0   ###### May cause problems? Check.
  if (!is.null(out_file)) {
    write.table(x,
                out_file,
                quote = FALSE,
                row.names = FALSE,
                sep = '\t')
  }
  
  # cycle through the condition pairs in the file and plot each pair
  for (i in seq_len(dim(repplot)[1])) {
    if(verbose) cat("--- PLOTTING REPLICATE PLOT ", i, ": ")
    
    # check if the replicate combination exists in the plots
    rep1_1 <-
      paste(repplot$condition1[i], repplot$rep1_1[i], sep = "_")
    rep1_2 <-
      paste(repplot$condition1[i], repplot$rep1_2[i], sep = "_")
    rep2_1 <-
      paste(repplot$condition2[i], repplot$rep2_1[i], sep = "_")
    rep2_2 <-
      paste(repplot$condition2[i], repplot$rep2_2[i], sep = "_")
    reps <- c(rep1_1, rep1_2, rep2_1, rep2_2)
    if (!any(!(reps %in% names(x)))) {
      # prep 1st replicate comparison for plot
      rep1 <-
        log2(x[, paste(repplot$condition1[i], 
                       repplot$rep1_1[i], sep = "_")] / x[, paste(repplot$condition2[i], 
                                                                  repplot$rep2_1[i], 
                                                                  sep = "_")])
      # prep 2nd replicate comparison for plot
      rep2 <-
        log2(x[, 
               paste(repplot$condition1[i], 
                     repplot$rep1_2[i], sep = "_")] / x[, paste(repplot$condition2[i], repplot$rep2_2[i], sep =
                                                                                           "_")])
      
      # remove NA pairs
      idx <-
        which(is.na(rep1) |
                is.na(rep2) | is.infinite(rep1) | is.infinite(rep2))
      rep1 <- rep1[-idx]
      rep2 <- rep2[-idx]
      
      if (length(rep1) > 1 & length(rep2) > 1) {
        reps.cor <-
          cor(rep1, rep2, use = "pairwise.complete.obs", method = "pearson")
        # set up a square plot centered at 0
        x.lim <- ceiling(max(abs(c(rep1, rep2)), na.rm = TRUE))
        y.lim <- c(-x.lim, x.lim)
        x.lim <- c(-x.lim, x.lim)
        
        # name axes labels
        y.label <-
          paste0(
            repplot$condition1[i],
            " vs. ",
            repplot$condition2[i],
            "  (",
            repplot$rep1_1[i],
            "/",
            repplot$rep2_1[i],
            ")"
          )
        x.label <-
          paste0(
            repplot$condition1[i],
            " vs. ",
            repplot$condition2[i],
            "  (",
            repplot$rep1_2[i],
            "/",
            repplot$rep2_2[i],
            ")"
          )
        # make plot name
        plot.name <-
          paste(
            repplot$condition1[i],
            " vs "  ,
            repplot$condition2[i],
            " R = ",
            round(reps.cor, 3),
            sep = ""
          )
        plot.name2 <-
          paste(
            repplot$condition1[i],
            " vs "  ,
            repplot$condition2[i],
            " R ",
            round(reps.cor, 3),
            sep = ""
          )
        tmp <- data.frame(rep1, rep2, stringsAsFactors = FALSE)
        p <- ggplot(tmp, aes(x = rep1, y = rep2)) +
          geom_point() +
          xlim(x.lim[1], x.lim[2]) +
          ylim(x.lim[1], x.lim[2]) +
          ggtitle(plot.name) +
          labs(x = x.label, y = y.label)
        
        if (!is.null(out_file)) {
          pdf_nameout <-
            paste(
              dirname(out_file),
              "/",
              gsub(" ", "_", plot.name2) ,
              "_",
              repplot$rep1_1[i],
              "_",
              repplot$rep1_2[i],
              "-",
              prot_exp,
              ".pdf",
              sep = ""
            )
          ggsave(
            filename = pdf_nameout,
            plot = p,
            width = 10,
            height = 10
          )
          
          if(verbose) cat(pdf_nameout, "\n")
        }
      } else{
        if(verbose) cat(
          "WARNING: not enough data for correlation analysis\n"
        )
      }
      } else{
        warning(
          "--- REPLICATE PLOT ",
          i,
          " NOT MADE -- MISSING DATA FROM ",
          paste(" ", reps[!(reps %in% names(x))], "\n", collapse = "")
        )
    }
  }
  }
