#' @rawNamespace import(biomaRt, except = c(select, getSequence))
#' @import bit64
#' @import corrplot
#' @import data.table
#' @import getopt
#' @import ggdendro
#' @import ggplot2
#' @import grid
#' @import limma
#' @import MSstats
#' @import org.Hs.eg.db
#' @import pheatmap
#' @rawNamespace import(plotly, except = c(last_plot, mutate, arrange, rename, summarise))
#' @import plyr
#' @import RColorBrewer
#' @import reshape2
#' @rawNamespace import(seqinr, except = c(zscore, count, a))
#' @import shiny
#' @import stringr
#' @import VennDiagram
#' @import yaml
#' @importFrom graphics pairs plot
#' @importFrom grDevices colorRampPalette dev.off pdf
#' @importFrom stats aggregate as.dendrogram cor dist fisher.test hclust kmeans median order.dendrogram
#' @importFrom utils combn read.delim write.table


# ------------------------------------------------------------------------------
#' @title Main Function
#' @description Main function running all the selected options
#' @param yaml_config_file The yaml file name and location
#' @return All the selected options
#' @keywords main, driver, function
#' artms_main()
#' @export
artms_main <- function(yaml_config_file){
  cat(">> RUNNING artMS. LOADING CONFIGURATION FILE...\n")
  config <- yaml.load_file(yaml_config_file)
  # process MaxQuant data, link with keys, and convert for MSStats format
  if(config$data$enabled){
    cat(">> LOADING DATA\n")
    ## Found more bugs in fread (issue submitted to data.table on github by
    ## JVD but it was closed with the excuse that 'is was not reproducible'
    ## although he provided examples) 
    ## Not worth the compromise in data integrity just to save time 
    ## reading in data
    
    # CHECKING FOR SILAC EXPERIMENT
    if(!is.null(config$silac$enabled)){
      if(config$silac$enabled){
        output <- gsub(".txt","-silac.txt",config$files$data)
        data <- artms_SILACtoLong(config$files$data, output)
      }else{
        data <- read.delim(config$files$data, stringsAsFactors=F, sep='\t')
      }
    }else{
      data <- read.delim(config$files$data, stringsAsFactors=F, sep='\t')
    }
    
    data <- data.table(data)
    setnames(data, colnames(data), gsub('\\s','.',colnames(data)))
    
    keys <- read.delim(config$files$keys, stringsAsFactors=F, sep='\t')
    keys <- data.table(keys)
    
    if( !any(grepl("RawFile", names(data))) ){
      tryCatch(setnames(data, 'Raw.file', 'RawFile'), error=function(e) cat('Raw.file not found in the evidence file\n'))
    }
    if( !any(grepl("RawFile", names(keys))) ){
      tryCatch(setnames(keys, 'Raw.file', 'RawFile'), error=function(e) cat('Raw.file not found in the KEYS file (and it should crash)\n'))
    }
    
    cat('\tVERIFYING DATA AND KEYS\n')
    
    if(!'IsotopeLabelType' %in% colnames(data)){
      cat("------- + IsotopeLabelType not detected in evidence file! 
          It will be assumed that this is a label-free experiment 
          (adding IsotopeLabelType column with L value\n")
      data[,IsotopeLabelType:='L']
    }
    
    # HACK FOR SILAC DATA
    if(!is.null(config$silac$enabled)){
      if(config$silac$enabled){
        data$RawFile = paste(data$RawFile, data$IsotopeLabelType, sep='')
        keys$RawFile = paste(keys$RawFile, keys$IsotopeLabelType, sep='')
        keys$Run = paste(keys$IsotopeLabelType,keys$Run , sep='')
        data$IsotopeLabelType = 'L'
        keys$IsotopeLabelType = 'L'
        data <- artms_mergeMaxQDataWithKeys(data, keys, by = c('RawFile','IsotopeLabelType'))
      }else{
        data <- artms_mergeMaxQDataWithKeys(data, keys, by = c('RawFile','IsotopeLabelType'))
      }
    }else{
      data <- artms_mergeMaxQDataWithKeys(data, keys, by = c('RawFile','IsotopeLabelType'))
    }
    
    ## fix for weird converted values from fread
    data[Intensity<1,]$Intensity=NA 
    
    ## FILTERING : handles Protein Groups and Modifications
    if(config$filters$enabled) data_f = artms_filterData(data, config) else data_f=data
    
    ## FORMATTING IN WIDE FORMAT TO CREATE HEATMAPS
    if(!is.null(config$files$sequence_type)){
      cat(">> OLD CONFIGUATION FILE DETECTED : sequence_type DETECTED. 
          WARNING: RECOMMENDED TO ALWAYS USED modified HERE\n")
      if(config$files$sequence_type == 'modified') castFun = artms_castMaxQToWidePTM else castFun = artms_castMaxQToWide
      data_w = castFun(data_f)
    }else{
      data_w = artms_castMaxQToWidePTM(data_f)
    }
    
    ## HEATMAPS
    if(!is.null(config$files$sample_plots) && config$files$sample_plots){
      keys_in_data = keys[keys$RawFile %in% unique(data$RawFile),]
      artms_sampleCorrelationHeatmap(data_w = data_w, keys = keys_in_data, config = config) 
      artms_samplePeptideBarplot(data_f, config)
    }
  }
  
  ## MSSTATS
  if(config$msstats$enabled){
    
    if(is.null(config$msstats$msstats_input)){
      # Go through the old yaml version. 
      # Before "fractions" it was called "aggregation" in the config.yaml file
      if(!is.null(config$aggregation$enabled)){
        config$fractions$enabled <- config$aggregation$enabled
        config$fractions$aggregate_fun <- config$aggregation$aggregate_fun
      }
      
      # DEPRECATED OPTION: in older versions the type of sequence 
      # could be selected (either modified or unmodified). 
      if(is.null(config$files$sequence_type)){
        config$files$sequence_type <- 'modified'
      }
      
      dmss <- artms_getMSstatsFormat(data_f, config$fractions$enabled, config$files$data, config$fractions$aggregate_fun)
      
      ## DEPRECATED : Make sure there are no doubles !!
      ## doubles could arise when protein groups are being kept and the same 
      ## peptide is assigned to a unique Protein. Not sure how this is possible 
      ## but it seems to be like this in maxquant output. 
      ## A possible explanation is that these peptides have different 
      ## retention times (needs to be looked into)
      ## dmss <- data.frame(dmss[,j=list(ProteinName=paste(ProteinName,collapse=';'),Intensity=median(Intensity, na.rm=T)),by=c('PeptideSequence','ProductCharge','PrecursorCharge','FragmentIon','IsotopeLabelType','Run','BioReplicate','Condition')])
      
    }else{
      cat(sprintf("\tREADING PREPROCESSED\t%s\n", config$msstats$msstats_input)) 
      dmss <- read.delim(config$msstats$msstats_input, stringsAsFactors=F, sep='\t')
      dmss <- data.table(dmss)
    }
    
    # Read in contrasts file
    # contrasts = read.delim(config$files$contrasts, stringsAsFactors=F)
    contrasts <- artms_writeContrast(config$files$contrasts)
    # make sure the column names are in alphabetical order before continuing
    contrasts = as.matrix( contrasts[,order(dimnames(contrasts)[[2]], decreasing=F)] )
    results = artms_runMSstats(dmss, contrasts, config)
  }
  
  ## ANNOTATING RESULT FILE
  if(config$output_extras$enabled){
    if(!config$msstats$enabled) results = read.delim(config$output_extras$msstats_output, stringsAsFactors=F)
    artms_writeExtras(results$ComparisonResult, config)
  }
  
  cat(">> ANALYSIS COMPLETE! HAVE A NICE DAY :)\n")
}


