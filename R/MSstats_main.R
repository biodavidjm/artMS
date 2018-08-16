# annotate the files based on the Uniprot accession id's
Extras.annotate <- function(results, output_file=opt$output, uniprot_ac_col='Protein', group_sep=';', uniprot_dir = '~/github/kroganlab/source/db/', species='HUMAN'){
  cat(">> ANNOTATING\n")
  
  # remove unnamed proteins that are listed as ""
  if(length(which(results$Protein==""))>0)  results <- results[-which(results$Protein==""),]
  
  # read in all the annotation files from the uniprot_db directory
  species_split = unlist(strsplit(species, "-"))
  Uniprot = NULL
  for(org in species_split){
    cat(sprintf("\tLOADING %s\n",org))
    tmp = read.delim(sprintf("%s/uniprot_protein_descriptions_%s.txt",uniprot_dir,org), stringsAsFactors=F, quote="")    
    if(is.null(Uniprot)){
      Uniprot = as.data.frame(tmp)
    }else{
      Uniprot = rbind(Uniprot, tmp)  
    }
  }
  
  # get list of all unique prey entries in this file. Keep 'group_sep' in mind.
  preys <- unique(results$Protein)
  preys <- preys.original <- data.frame(prey = preys, idx = 1:length(preys), stringsAsFactors=F)
  # split apart all the preys and index them so we can piece them back together when merging
  preys <- do.call(rbind, apply(preys, 1, function(y){ data.frame(prey = unlist(strsplit(y[1], ";")), idx = as.numeric(y[2]), stringsAsFactors = F) } ))
  
  # annotate all the preys wiht the Uniprot info
  preys <- merge( preys, Uniprot[,c("Entry","Entry.name","Protein.names","Gene.names")], by.x ="prey", by.y="Entry", all.x=T)
  # aggregate all the preys on the indexes so we can merge with the original data
  
  # merge protein name
  tmp <- aggregate(data = preys[,c("prey",'idx','Entry.name')], .~idx, paste, collapse=";")
  names(tmp) <- c('idx','uniprot_ac','Protein_name')
  preys.new <- merge(preys.original, tmp, by='idx', all.x=T)
  
  # merge protein description
  tmp <- aggregate(data = preys[,c('idx','Protein.names')], .~idx, paste, collapse=";")
  names(tmp) <- c('idx', 'Protein_desc')
  preys.new <- merge(preys.new, tmp, by='idx', all.x=T)
  
  # merge protein description
  preys$Gene.names <- gsub(" .*","", preys$Gene.names)
  tmp <- aggregate(data = preys[,c('idx','Gene.names')], .~idx, paste, collapse=";")
  names(tmp) <- c('idx','Gene.names')
  preys.new <- merge(preys.new, tmp, by='idx',all.x=T)
  
  # merge the annotations all back into the original data
  results_out <- merge(results, preys.new, by.x="Protein", by.y="prey", all.x=T)
  results_out$idx = c()
  
  # alert user of any unmapped proteins
  unmapped = unique(results_out[is.na(results_out$uniprot_ac),"Protein"]) 
  cat('UNMAPPED PROTEINS\t', length(unmapped), '\n')
  cat('\t',paste(unmapped,collapse='\n\t'),'\n')
  write.table(results_out, file=output_file, sep='\t', quote=F, row.names=F, col.names=T)
  cat(">> ANNOTATING COMPLETE!\n")
  return(results_out)
}

# ------------------------------------------------------------------------------
#' @title Write extras
#' @description Extras after MSstats, as annotations, volcano plots, heatmaps
#' @param results MSstats results
#' @param config The configuration object (yaml)
#' @return Extras as selected in the yaml file, including:
#' - volcano plot (pdf)
#' - Adding annotations (gene symbol based on uniprot)
#' - 
#' @keywords extras, annotations, volcano
#' writeExtras()
#' @export
writeExtras <- function(results, config){
  
  if(length(results)==0 | !exists('results')){
    stop("ERROR!! NO RESULTS FOUND TO ANNOTATE!")
  }
  
  # Annotation 
  if(config$output_extras$annotate & is.null(config$filters$modifications) ){
    results_ann <- Extras.annotate(results, output_file=config$files$output, uniprot_ac_col='Protein', group_sep=';', uniprot_dir = config$output_extras$annotation_dir, species=config$output_extras$species)
  }else{
    if( !is.null(config$filters$modifications) ) cat("\tSITES NEED TO BE MAPPED BACK TO PROTEINS BEFORE ANNOTATING.\n")
    results_ann = results
    if( !is.null(config$output_extras$msstats_output)){
      config$files$output = config$output_extras$msstats_output
    }else{
      cat("\tNO PREVIOUS MSSTAT OUTPUT FILE NOTED. USING CURRENT RESULTS FOR EXTRAS.\n")
    }
  }
  
  lfc_lower = as.numeric(unlist(strsplit(config$output_extras$LFC,split=" "))[1])
  lfc_upper = as.numeric(unlist(strsplit(config$output_extras$LFC,split=" "))[2])
  ## select subset of labels for heatmap and volcan plots
  selected_labels = config$output_extras$comparisons
  if(is.null(selected_labels) || selected_labels=='all') selected_labels='*'
  
  # remove the Inf, -Inf log2FC hits. 
  results_ann <- results_ann[!is.infinite(results_ann$log2FC),]
  
  ## select data points  by LFC & FDR criterium in single condition and adding corresponding data points from the other conditions
  sign_hits = significantHits(results_ann,labels=selected_labels,LFC=c(lfc_lower,lfc_upper),FDR=config$output_extras$FDR)
  if( dim(sign_hits)[1] == 0 ) stop("NO SIGNIFICANT HITS DETECTED IN THIS EXPERIMENT. ABORTING PLOTS.\n")
  sign_labels = unique(sign_hits$Label)
  cat(sprintf("\tSELECTED HITS FOR PLOTS WITH LFC BETWEEN %s AND %s AT %s FDR:\t%s\n",lfc_lower, lfc_upper, config$output_extras$FDR, nrow(sign_hits)/length(sign_labels))) 
  
  cat(paste("output_file: ", config$files$output, "\n"))
  
  ## REPRESENTING RESULTS AS HEATMAP
  if(config$output_extras$heatmap){
    ## plot heat map for all contrasts
    cat(">>   PLOTTING HEATMAP FOR ALL CONTRASTS\n")
    heat_labels = prettyPrintHeatmapLabels(uniprot_acs=sign_hits$Protein,uniprot_ids=sign_hits$name, gene_names=sign_hits$Gene.names)
    heat_data_w = plotHeat(mss_F = sign_hits, out_file =  gsub('.txt','-sign.pdf',config$files$output), names=heat_labels, cluster_cols=config$output_extras$heatmap_cluster_cols, display = config$output_extras$heatmap_display)  
  }
  
  if(config$output_extras$volcano){
    cat(">>   PLOTTING VOLCANO PLOT\n")
    file_name = gsub('.txt','-volcano.pdf',config$files$output)
    volcanoPlot(results_ann[grep(selected_labels,results_ann$Label),], lfc_upper, lfc_lower, FDR=config$output_extras$FDR, file_name=file_name)  
  }
}


#' @title Main Function
#' @description Main function running all the selected options
#' @param opt the object with all the options
#' @return All the selected options
#' @keywords main, driver, function
#' main()
#' @export
main <- function(opt){
  cat(">> MSSTATS PIPELINE\n")
  config = tryCatch(yaml.load_file(opt$config), error = function(e) {cat(opt$config);break} )
  
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
        data <- MQutil.SILACToLong(config$files$data, output)
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
        data <- mergeMaxQDataWithKeys(data, keys, by = c('RawFile','IsotopeLabelType'))
      }else{
        data <- mergeMaxQDataWithKeys(data, keys, by = c('RawFile','IsotopeLabelType'))
      }
    }else{
      data <- mergeMaxQDataWithKeys(data, keys, by = c('RawFile','IsotopeLabelType'))
    }
    
    ## fix for weird converted values from fread
    data[Intensity<1,]$Intensity=NA 
    
    ## FILTERING : handles Protein Groups and Modifications
    if(config$filters$enabled) data_f = filterData(data, config) else data_f=data
    
    ## FORMATTING IN WIDE FORMAT TO CREATE HEATMAPS
    if(!is.null(config$files$sequence_type)){
      cat(">> OLD CONFIGUATION FILE DETECTED : sequence_type DETECTED. 
          WARNING: RECOMMENDED TO ALWAYS USED modified HERE\n")
      if(config$files$sequence_type == 'modified') castFun = castMaxQToWidePTM else castFun = castMaxQToWide
      data_w = castFun(data_f)
    }else{
      data_w = castMaxQToWidePTM(data_f)
    }
    
    ## HEATMAPS
    if(!is.null(config$files$sample_plots) && config$files$sample_plots){
      keys_in_data = keys[keys$RawFile %in% unique(data$RawFile),]
      sampleCorrelationHeatmap(data_w = data_w, keys = keys_in_data, config = config) 
      samplePeptideBarplot(data_f, config)
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
      
      dmss <- getMSstatsFormat(data_f, config$fractions$enabled, config$files$data, config$fractions$aggregate_fun)
      
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
    contrasts <- writeContrast(config$files$contrasts)
    # make sure the column names are in alphabetical order before continuing
    contrasts = as.matrix( contrasts[,order(dimnames(contrasts)[[2]], decreasing=F)] )
    results = runMSstats(dmss, contrasts, config)
  }
  
  ## ANNOTATING RESULT FILE
  if(config$output_extras$enabled){
    if(!config$msstats$enabled) results = read.delim(config$output_extras$msstats_output, stringsAsFactors=F)
    writeExtras(results$ComparisonResult, config)
  }
  
  cat(">> ANALYSIS COMPLETE! HAVE A NICE DAY :)\n")
}


## CONFIG LOADING

# opt = getopt(spec = spec, opt = commandArgs(TRUE), command = get_Rscript_filename(), usage = FALSE, debug = FALSE)
# opt = c(opt, config='~/experiments/artms/technical_replicas/configTR.yaml')
