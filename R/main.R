#' @rawNamespace import(AnnotationDbi, except = c(head))
#' @rawNamespace import(biomaRt, except = c(select, getSequence)) #delete?
#' @import bit64
#' @import circlize
#' @importFrom cluster pam
#' @import ComplexHeatmap
#' @importFrom corrplot corrplot
#' @importFrom dplyr mutate desc count arrange desc
#' @rawNamespace import(data.table, except = c(melt))
#' @importFrom factoextra fviz_pca_var fviz_contrib fviz_pca_ind 
#' fviz_nbclust get_dist fviz_silhouette fviz_cluster
#' @importFrom FactoMineR PCA
#' @import getopt
#' @import ggdendro
#' @import ggplot2
#' @importFrom gplots heatmap.2
#' @import ggrepel
#' @import gProfileR
#' @importFrom graphics pairs plot barplot hist lines par panel.smooth rect strwidth text
#' @importFrom grDevices colorRampPalette dev.off pdf
#' @import grid
#' @import limma
#' @import MSstats
#' @import openxlsx
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @importFrom PerformanceAnalytics chart.Correlation
#' @import pheatmap
#' @rawNamespace import(plotly, except = c(last_plot, mutate, arrange, 
#' rename, summarise, select, add_heatmap))
#' @importFrom plyr ddply summarise
#' @import RColorBrewer
#' @importFrom reshape2 melt
#' @rawNamespace import(seqinr, except = c(zscore, count, a))
#' @importFrom stats aggregate as.dendrogram cor dist fisher.test hclust prcomp quantile sd
#' kmeans median order.dendrogram phyper as.dist complete.cases
#' @import stringr
#' @importFrom tidyr unnest
#' @import UpSetR
#' @importFrom utils combn read.delim sessionInfo write.table setTxtProgressBar 
#' txtProgressBar head globalVariables
#' @import VennDiagram
#' @import yaml

utils::globalVariables(
  c(
    "...",
    "..count..",
    "AbMean",
    "Abundance",
    "ABUNDANCE",
    "adj.pvalue",
    "artms_data_corum_mito_database",
    "artms_data_pathogen_LPN",
    "artms_data_pathogen_TB",
    "artms_config",
    "bin.all",
    "bin.condition",
    "BioReplicate",
    "bioreplicate",
    "category",
    "charge",
    "cluster",
    "Comparison",
    "ComplexName",
    "condition",
    "Condition",
    "Contaminant",
    "ENTREZID",
    "evidencekeys",
    "experiment",
    "Experiment",
    "FEATURE",
    "fraction",
    "fxDx",
    "FxOverSamp",
    "Gene",
    "GENENAME",
    "GROUP_ORIGINAL",
    "iLog2FC",
    "imputedDFext",
    "installed.packages",
    "IntDetection",
    "intensity",
    "Intensity",
    "Ions",
    "Ions.mean",
    "Ions.sem",
    "iPvalue",
    "isotope.patterns",
    "isotope.patterns.sequenced..z.1.",
    "IsotopeLabelType",
    "isPtm",
    "Label",
    "log2FC",
    "log2fc_file",
    "m.z",
    "mc",
    "missed.cleavages",
    "MODIFICATION",
    "Modifications",
    "modified.sequence",
    "ms",
    "ms.ms",
    "ms.ms.count",
    "ms.ms.identified....",
    "MSMS.counts",
    "N",
    "nsamples",
    "num.IsotopePatterns.mean",
    "num.IsotopePatterns.sem",
    "num.IsotopePatternsSeq.mean",
    "num.IsotopePatternsSeq.sem",
    "num.MS1.mean",
    "num.MS1.sem",
    "num.MS2.mean",
    "num.MS2.sem",
    "output_dir",
    "oxidation..m.",
    "pathogen.ids",
    "PC1",
    "PC2",
    "pct.MS2Id.mean",
    "pct.MS2Id.sem",
    "pct.OxM",
    "pCV",
    "pDX",
    "PeakWidth.mean",
    "pearson",
    "Peptides",
    "Peptides.mean",
    "Peptides.sem",
    "potential.contaminant",
    "Prey",
    "prot_names",
    "Protein",
    "PROTEIN",
    "proteins",
    "Proteins",
    "Proteins.mean",
    "Proteins.sem",
    "PSMs",
    "PSMs.mean",
    "PSMs.sem",
    "ptm_site",
    "PTMone",
    "PTMsite",
    "pvalue",
    "quantile",
    "RawFile",
    "rect",
    "BiorepCount",
    "CondCount",
    "res_index",
    "retention.length",
    "reverse",
    "sample_name",
    "Species",
    "species",
    "strwidth",
    "SUBJECT_ORIGINAL",
    "sumInt",
    "SYMBOL",
    "text",
    "TR",
    "tr1",
    "tr2",
    "type",
    "uncalibrated.mass.error..ppm.",
    "uniprot_ac",
    "uniprot_id",
    "value",
    "variable",
    "ymax",
    "ymin"
  )
)

# ------------------------------------------------------------------------------
#' @title Relative quantification using MSstats
#'
#' @description Relative quantification using MSstats including:
#' - plots
#' - quantifications (log2fc, pvalues, etc)
#' - normalized abundance values
#' @param yaml_config_file (char) The yaml file name and location
#' @param verbose (logical) `TRUE` (default) shows function messages
#' @return The relative quantification of the conditions and comparisons
#' specified in the keys/contrast file resulting from running MSstats, in
#' addition to quality control plots (if selected)
#' @keywords main, driver, function
#' @examples \donttest{
#' artmsQuantification("artms-ab-config.yaml")
#' }
#' @export
artmsQuantification <- function(yaml_config_file,
                                 verbose = TRUE) {
  
  if(verbose){
    message(" WELCOME to artMS (Analytical R Tools for Mass Spectrometry) ")
    message("============================================================  ")
    message(">> LOADING CONFIGURATION FILE... ")
  }
  
  config <- yaml.load_file(yaml_config_file)
  
  # CHECK POINT: DO THE FILES EXIST?
  if(!file.exists(config$files$contrasts)){
    stop("the file ", config$files$contrasts, " does not exist! ")
  }
  
  if(!file.exists(config$files$evidence)){
    stop("the file ", config$files$evidence, " does not exist! ")
  }
  
  if(!file.exists(config$files$keys)){
    stop("the file ", config$files$keys, " does not exist! ")
  }
  
  if(!(grepl("\\.txt$", config$files$output))){
    stop("the file ", config$files$output, " must have extension .txt " )
  }
  
  
  # LET'S HELP THE DISTRACTED USER
  if (!(is.null(config$data$filters$modification))) {
    config$data$filters$modification <-
      toupper(config$data$filters$modification)
  }
  
  # Quality Control
  if (config$qc$basic) {
    artmsQualityControlEvidenceBasic(
      evidence_file = config$files$evidence,
      keys_file = config$files$keys,
      prot_exp = toupper(config$data$filters$modifications),
      fractions = config$data$fractions$enabled)
  }
  
  if (config$qc$extended) {
    artmsQualityControlEvidenceExtended  (
      evidence_file = config$files$evidence,
      keys_file = config$files$keys)
  }
  
  # process MaxQuant data, link with keys, and convert for MSStats format
  if (config$data$enabled) {
    if(verbose)
      message(" QUANTIFICATION: LOADING DATA----------------- ")
    ## Found more bugs in fread (issue submitted to data.table on github by
    ## JVD but it was closed with the excuse that 'is was not reproducible'
    ## although he provided examples)
    ## Not worth the compromise in data integrity just to save time
    ## reading in data
    
    # CHECKING FOR SILAC EXPERIMENT
    if (!is.null(config$data$silac$enabled)) {
      if (config$data$silac$enabled) {
        output <- gsub(".txt", "-silac.txt", config$files$evidence)
        x <- artmsSILACtoLong(config$files$evidence,
                                  output,
                                  verbose = verbose)
      } else{
        x <- .artms_checkIfFile(config$files$evidence)
        x <- .artms_checkRawFileColumnName(x)
      }
    } else{
      x <- .artms_checkIfFile(config$files$evidence)
      x <- .artms_checkRawFileColumnName(x)
    }
    
    x <- data.table(x)
    
    keys <- .artms_checkIfFile(config$files$keys)
    keys <- .artms_checkRawFileColumnName(keys)
    
    keys <- data.table(keys)
    
    # Let's make sure that the contrast file is right
    if (config$msstats$enabled) {
      # Read in contrast file
      contrasts <-
        .artms_writeContrast(
          contrast_file = config$files$contrasts, 
          all_conditions = unique(as.character(keys$Condition)))
    }
    if(verbose) message('\tVERIFYING DATA AND KEYS ')
    
    if (!'IsotopeLabelType' %in% colnames(x)) {
      if(verbose) message(
        "------- + IsotopeLabelType not detected in evidence file!
        It will be assumed that this is a label-free experiment
        (adding IsotopeLabelType column with L value) "
      )
      x[, IsotopeLabelType := 'L']
    }
    
    # HACK FOR SILAC DATA
    if (!is.null(config$data$silac$enabled)) {
      if (config$data$silac$enabled) {
        x$RawFile = paste(x$RawFile, x$IsotopeLabelType, sep = '')
        keys$RawFile = paste(keys$RawFile, keys$IsotopeLabelType, sep =
                               '')
        keys$Run = paste(keys$IsotopeLabelType, keys$Run , sep = '')
        x$IsotopeLabelType = 'L'
        keys$IsotopeLabelType = 'L'
        x <-
          artmsMergeEvidenceAndKeys(x, 
                                     keys, 
                                     by = c('RawFile', 'IsotopeLabelType'),
                                     verbose = verbose)
      } else{
        x <-
          artmsMergeEvidenceAndKeys(x, 
                                     keys, 
                                     by = c('RawFile', 'IsotopeLabelType'),
                                     verbose = verbose)
      }
    } else{
      x <-
        artmsMergeEvidenceAndKeys(x, 
                                   keys, 
                                   by = c('RawFile', 'IsotopeLabelType'),
                                   verbose = verbose)
    }
    
    ## fix for weird converted values from fread
    x[Intensity < 1, ]$Intensity <- NA
    
    ## FILTERING : handles Protein Groups and Modifications
    if (config$data$filters$enabled){
      data_f <- .artms_filterData(x = x, 
                                  config = config, 
                                  verbose = verbose)
    }else{
      data_f <- x
    }
    
    ## FORMATTING IN WIDE FORMAT TO CREATE HEATMAPS
    if (!is.null(config$files$sequence_type)) {
      if(verbose)
      message(">> OLD CONFIGUATION FILE DETECTED : sequence_type DETECTED.
        WARNING: RECOMMENDED TO ALWAYS USED modified HERE ")
      if (config$files$sequence_type == 'modified'){
        castFun = .artms_castMaxQToWidePTM
      }else{
        castFun = .artms_castMaxQToWide
      }
      data_w = castFun(data_f)
    } else{
      data_w = .artms_castMaxQToWidePTM(data_f)
    }
    
    ## HEATMAPS
    if (!is.null(config$data$sample_plots) &&
        config$data$sample_plots) {
      keys_in_data <- keys[keys$RawFile %in% unique(x$RawFile), ]
      .artms_sampleCorrelationHeatmap(data_w = data_w,
                                      keys = keys_in_data,
                                      config = config)
      .artms_samplePeptideBarplot(data_f, config)
    }
  }
  
  ## MSSTATS
  if (config$msstats$enabled) {

    # Read in contrast file
    contrasts <-
      .artms_writeContrast(
        contrast_file = config$files$contrasts, 
        all_conditions= unique(as.character(keys$Condition)))
    
    selectedConditions <- as.character(colnames(contrasts))
    
    if (is.null(config$msstats$msstats_input)) {
      # Check point to prevent MSstats crashed in the number of conditions
      # in the comparisons is not the same than the one in the keys file
      data_f <- data_f[which(data_f$Condition %in% selectedConditions),]
      dmss <- .artms_getMSstatsFormat(data_f = data_f,
                                      fraction = config$data$fractions$enabled,
                                      output_name = config$files$evidence,
                                      funfunc = "sum")
    } else {
      if(verbose) message(sprintf("\tREADING PREPROCESSED\t%s ",
        config$msstats$msstats_input))
      dmss <- read.delim(config$msstats$msstats_input,
                         stringsAsFactors = FALSE,
                         sep = '\t')
      dmss <- data.table(dmss)
    }
    results <- .artms_runMSstats(dmss, contrasts, config,
                                 verbose = verbose)
  }
  
  ## ANNOTATING RESULT FILE
  if (config$output_extras$enabled) {
    if (!config$msstats$enabled){
      message("MSstats was not enabled, therefore <output_extras> cannot be done!")
    }else{
      .artms_writeExtras(results$ComparisonResult, config)
    }
  }
  if(verbose) message(" ANALYSIS COMPLETED ")
}

# ------------------------------------------------------------------------------
#' @title Write out a template file of the artMS configuration file (yaml)
#' 
#' @description Creates a template file of the artMS configuration file, which
#' is required to run `artmsQuantification`. Check `?artms_config` and the 
#' vignettes to find out more about the details of the structure of the file
#' and how to fill it up
#' @param config_file_name (char) The name for the configuration file. It must
#' have a `.yaml` extension. If `NULL`, it returns the config as a yaml object
#' @param verbose (logical) `TRUE` (default) shows function messages
#' @return A file (or yaml data object) of the artMS configuration file
#' @keywords config, yaml
#' @examples 
#' config_empty <- artmsWriteConfigYamlFile(config_file_name = NULL)
#' @export
artmsWriteConfigYamlFile <- function(
  config_file_name = "artms_config_file.yaml",
  verbose = TRUE){
  
  if(!is.null(config_file_name)){
    if(grepl("\\.yaml", config_file_name)){
      write_yaml(x = artms_config, file = config_file_name )
      if(verbose) message(">> File ",config_file_name,"is out ")
    }else{
      stop("The <config_file_name> must have the extension .yaml")
    }
  }else{
    return(artms_config)
  }

}