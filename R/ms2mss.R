#' @import data.table

# Format MaxQuant file to MSStats format ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Convert MaxQuant evidence file to MSStats compatible format
#' @description This will take a MaxQuant searched evidence data and convert it to a format compatible with MSStats.
#' @param evidence_file The filepath to the MaxQuant searched data (evidence) file (txt tab delimited file).
#' @param keys_file The filepath to the keys file used with MSStats (txt tab delimited file).
#' @param protein_groups How to handle protein groups. Default = 'ignore'
#'   \itemize{
#'     \item{'remove'}{: This will remove the protein groups entirely from the data set.}
#'     \item{'explode'}{: This will create a new row for each protein in the group, copying the MS data for all proteins.}
#'     \item{'ignore'}{: This will treate the protein group as if it were a normal protein unique identifier.}
#' }
#' @param contaminants Whether or not to filter out MaxQuant labeled contaminants. Default = TRUE
#' @param modification What kind of modified sites are in the data. Options : PH, UB. Default = NULL
#' @param plot.heatmap Whether or not to output a hierarchical clustered heatmap of the samples. Default = FALSE
#' @param out_file The filepath to where you want the heatmap to be output to as a pdf.
#' @param aggregate_intensities Whether or not to aggregate the intensities based on the unique ppairing of RawFile, Proteins, Sequence, Charge, and IsotopeLabelType. Default = FALSE
#' @param aggregate_fun The aggregate to apply to the Intensities. Default is 'max'.
#' @keywords MaxQuant, evidence
#' mq2mss()
#' @export
mq2mss <- function(evidence_file, keys_file, protein_groups = "ignore", contaminants = TRUE, modification = NULL, plot.heatmap = FALSE, 
    output_file = NULL, aggregate_intensities = FALSE, aggregate_fun = NULL) {
    
    cat(">> LOADING DATA\n")
    # read in evidence file
    data <- checkIfFile(evidence_file, is.evidence = T)
    data.table::setnames(data, colnames(data), gsub("\\s", ".", colnames(data)))
    keys = checkIfFile(keys_file)
    
    ## the following lines were added to integrate the Label with the Filename when using multiple labels (e.g. H/L) currently we
    ## leave this in because the MSstats discinction on labeltype doesn't work see ISSUES
    ## https://github.com/everschueren/RMSQ/issues/1
    
    tryCatch(data.table::setnames(data, "Raw.file", "RawFile"), error = function(e) cat("Raw.file not found\n"))
    tryCatch(data.table::setnames(keys, "Raw.file", "RawFile"), error = function(e) cat("Raw.file not found\n"))
    
    cat("\tVERIFYING DATA AND KEYS\n")
    if (!"IsotopeLabelType" %in% colnames(data)) 
        data[, `:=`(IsotopeLabelType, "L")]
    
    data = mergeMaxQDataWithKeys(data, keys, by = c("RawFile", "IsotopeLabelType"))
    data$RawFile = paste(data$RawFile, data$IsotopeLabelType, sep = "")
    keys$RawFile = paste(keys$RawFile, keys$IsotopeLabelType, sep = "")
    keys$Run = paste(keys$IsotopeLabelType, keys$Run, sep = "")
    data$IsotopeLabelType = "L"
    keys$IsotopeLabelType = "L"
    data[Intensity < 1, ]$Intensity = NA  ## fix for weird converted values from fread
    
    ## end hacks for SILAC
    
    ## FILTERING : handles MQ contaminants, Protein Groups, and Modifications
    config = c()
    config$filters$protein_groups = protein_groups
    config$filters$contaminants = contaminants
    config$filters$modification = modification
    config$files$data = output_file
    if (!is.null(config$filters$modification) | (config$filters$contaminants != F) | (config$filters$protein_groups != "ignore")) 
        data_f = filterData(data, config) else data_f = data  #!!!!!!!
    
    ## FORMATTING IN WIDE FORMAT FOR NORMALIZATION PURPOSES
    if (!is.null(config$filters$modification)) 
        castFun = castMaxQToWidePTM else castFun = castMaxQToWide  #!!!!!! SOMETHING WEIRD HERE
    data_w = castFun(data_f)
    
    # Create heatmap
    if (plot.heatmap) {
        keys_in_data = keys[keys$RawFile %in% unique(data$RawFile), ]
        if (is.null(output_file)) 
            stop("No output file designated for heatmap! Please enter a file path for the heatmap to be saved.\n")
        config$files$output = output_file
        sampleCorrelationHeatmap(data_w = data_w, keys = keys_in_data, config = config)
        samplePeptideBarplot(data_f, config)
    }
    
    ## AGGREGATION Aggregate Intensitis based on unique RawFileCombined + Proteins + Sequence + Charge + IsotopeLabelType
    ## combinations
    if (aggregate_intensities) {
        if (!is.null(aggregate_fun)) {
            config$aggregation$aggregate_fun = aggregate_fun
        } else {
            config$aggregation$aggregate_fun = NA
        }
        res = aggregateData(data_w, keys, config)
        data_w = res$data_w_agg
        keys = res$keys_agg
    }
    
    ## MSSTATS Convert to MSStats format
    dmss = data.table(convertDataLongToMss(data_w, keys, config))
    ## make sure there are no doubles !!  doubles could arise when protein groups are being kept and the same peptide is assigned to
    ## a unique Protein. Not sure how this is possible but it seems to be like this in maxquant output. A possible explanation is
    ## that these peptides have different retention times (needs to be looked into)
    dmss = data.frame(dmss[, j = list(ProteinName = paste(ProteinName, collapse = ";"), Intensity = median(Intensity, na.rm = T)), 
        by = c("PeptideSequence", "ProductCharge", "PrecursorCharge", "FragmentIon", "IsotopeLabelType", "Run", "BioReplicate", 
            "Condition")])
    return(dmss)
}




















