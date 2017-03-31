#! /usr/local/bin/RScript --vanilla
#' @import data.table
#'
###############################
## FILE AND LIB LOADING #######

# cat("\nYou are using RMSQ v3\n\n")
# cat(">> LOADING EXTERNAL FILES AND LIBRARIES\n")

# suppressMessages(library(getopt))
# suppressMessages(library(yaml))
# suppressMessages(library(biomaRt))

# # set source directory
# args <- commandArgs(trailingOnly = F)
# scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))

## load all external files
# source("MSstats_functions.R")


#########################
## MAIN FUNCTIONS #######

filterData = function(data, config){
  cat(">> FILTERING\n")
  if(config$filters$protein_groups == 'remove'){
    cat("\tPROTEIN GROUPS\tREMOVE\n")
    data_f = removeMaxQProteinGroups(data)
  }else if(config$filters$protein_groups == 'explode'){
    cat("\tPROTEIN GROUPS\tEXPLODE\n")
    data_f = explodeMaxQProteinGroups(data)
  }else{
    cat("\tPROTEIN GROUPS\tIGNORE\n")
    data_f = data
  }
  if(config$filters$contaminants){
    cat("\tCONTAMINANTS\tREMOVE\n")
    data_f = filterMaxqData(data_f)
  }
  if(!is.null(config$filters$modification)){
    cat(sprintf("\tMODIFICATIONS\t%s\n",config$filters$modification))
    if(config$filters$modification == 'UB'){
      data_f = data_f[Modifications %like% 'GlyGly']
    }else if(config$filters$modification == 'PH'){
      data_f = data_f[Modifications %like% 'Phospho']
    }
  }
  return(data_f)
}

aggregateData = function(data_w, keys, config, castFun){
  cat(">>AGGREGATING TECHNICAL REPEATS\n")
  ##  we need to convert to long format for more easy aggregation
  data_l = meltMaxQToLong(data_w, na.rm = T)
  keysagg = flattenKeysTechRepeats(keys)
  keysmapping = merge(keys, keysagg[,c('BioReplicate','RawFile')], by='BioReplicate')
  setnames(keysmapping,c('RawFile.x','RawFile.y'),c('RawFile','RawFileCombined'))
  data_l_combined = merge(data_l, keysmapping, by=c('RawFile','IsotopeLabelType'))
  aggregate_fun = tryCatch(eval(parse(text=config$aggregation$aggregate_fun)),
                           error = function(e) cat("\tWARNING: argument for aggregate_fun not a valid R function - defaulting to:\t 'max'"),
                           finally=max)

  data_l_combined_agg = data.table(aggregate(Intensity ~ RawFileCombined + Proteins + Sequence + Charge + IsotopeLabelType, FUN=sum, data=data_l_combined))
  setnames(data_l_combined_agg,'RawFileCombined','RawFile')
  data_w_agg = castMaxQToWide(data_l_combined_agg)
  return(list(data_w_agg = data_w_agg, keys_agg = keysagg))
}

## returns data tabel in wide format
normalizeData = function(data_w, config){
  cat(">> NORMALIZING\n")

  if(grepl('scale|quantile|cyclicloess',config$normalization$method)){
    cat(sprintf("\tNORMALIZATION\t%s\n",config$normalization$method))
    data_fn = normalizeSingle(data_w=data_w, NORMALIZATION_METHOD=config$normalization$method)
  }else if(grepl('reference',config$normalization$method) && !is.null(config$normalization$reference) && config$normalization$reference %in% unique(data_w$Proteins)){

    ## CURRENTLY UING MSSTATS METHODS FOR REFERENCE NORMALIZATION
    data_fn = data_w

#     cat(sprintf("\tNORMALIZATION\tTO REFERENCE\t%s\n",config$normalization$reference))
#     ref_files = keys[keys$NormalizationGroup == 'REFERENCE', 'RawFile']
#     data_l_ref = data_l[data_l$Raw.file %in% ref_files & data_l$Intensity > 0 & is.finite(data_l$Intensity), ]
#     data_l_nonref = data_l[!(data_l$Raw.file %in% ref_files) & data_l$Intensity > 0 & is.finite(data_l$Intensity), ]
#     data_l_ref_n = normalizeToReference(data_l_ref=data_l_ref, ref_protein = config$normalization$reference, output_file = config$files$output)
#     data_fn= rbind(data_l_ref_n, data_l_nonref)
#     data_fn = castMaxQToWide(data_fn)

  }else{
    data_fn = data_w
  }

  return(data_fn)
}

runMSstats = function(dmss, contrasts, config){
  # plot the data BEFORE normalization
  if(grepl('before', config$msstats$profilePlots)){
    mssquant = dataProcess(raw = dmss,
                           normalization=F,
                           betweenRunInterferenceScore = config$msstats$interference,
                           fillIncompleteRows = T,
                           summaryMethod = config$msstats$summaryMethod,
                           censoredInt = config$msstats$censoredInt,
                           cutoffCensored = config$msstats$cutoffCensored,
                           MBimpute = config$msstats$MBimpute,
                           featureSubset=config$msstats$feature_subset)
    dataProcessPlots(data=mssquant, type="ProfilePlot", featureName="Peptide", address=gsub('.txt','-before',config$files$output))
    dataProcessPlots(data=mssquant, type="QCPlot", address=gsub('.txt','-before',config$files$output))
  }

  # Normalization
  if(!is.null(config$msstats$normalization_reference) & config$msstats$normalization_method == 'globalStandards'){  # if globalStandars is selected, must have a reference protein(s)
    normalization_refs = unlist(lapply(strsplit(config$msstats$normalization_reference, split = ','), FUN=trim))
    #mssquant = dataProcess(dmss, normalization=config$msstats$normalization_method, nameStandards=normalization_refs , fillIncompleteRows=T)
    mssquant = dataProcess(raw = dmss,
                           normalization=config$msstats$normalization_method,
                           nameStandards=normalization_refs,
                           betweenRunInterferenceScore = config$msstats$interference,
                           fillIncompleteRows = T,
                           summaryMethod = config$msstats$summaryMethod,
                           censoredInt = config$msstats$censoredInt,
                           cutoffCensored = config$msstats$cutoffCensored,
                           MBimpute = config$msstats$MBimpute,
                           featureSubset=config$msstats$feature_subset)
  }else{
    cat(sprintf('>> NORMALIZATION\t%s\n',config$msstats$normalization_method))
    #mssquant = dataProcess(dmss, normalization=config$msstats$normalization_method , fillIncompleteRows = F, betweenRunInterferenceScore = F)
    mssquant = dataProcess(raw = dmss,
                           normalization=config$msstats$normalization_method,
                           betweenRunInterferenceScore = config$msstats$interference,
                           fillIncompleteRows = T,
                           summaryMethod = config$msstats$summaryMethod,
                           censoredInt = config$msstats$censoredInt,
                           cutoffCensored = config$msstats$cutoffCensored,
                           MBimpute = config$msstats$MBimpute,
                           featureSubset=config$msstats$feature_subset)
  }

  # plot the data AFTER normalization
  if(grepl('after', config$msstats$profilePlots)){
    dataProcessPlots(data=mssquant, type="ProfilePlot", featureName="Peptide", address=gsub('.txt','-after',config$files$output))
    dataProcessPlots(data=mssquant, type="QCPlot", address=gsub('.txt','-after',config$files$output))
  }

  if(!all(levels(mssquant$GROUP_ORIGINAL) == colnames(contrasts))){
    stop(sprintf('\tERROR IN CONTRAST COMPARISON: GROUP LEVELS DIFFERENT FROM CONTRASTS FILE\n\tGROUP LEVELS\t%s\n\tCONTRASTS FILE\t%s\n', paste(levels(mssquant$GROUP_ORIGINAL),collapse=','),paste(colnames(contrasts),collapse=',')))
  }

  # protein sample/group quantification
  write.table(quantification(mssquant), file=gsub('.txt','-mss-sampleQuant.txt',config$files$output), eol="\n", sep="\t", quote=F, row.names=F, col.names=T)

  cat(sprintf('\tFITTING CONTRASTS:\t%s\n',paste(rownames(contrasts),collapse=',')))
  write.table(mssquant$ProcessedData, file=gsub('.txt','-mss-normalized.txt',config$files$output), eol="\n", sep="\t", quote=F, row.names=F, col.names=T)
  #results = groupComparison(data = mssquant, contrast.matrix = contrasts, labeled = as.logical(config$msstats$labeled), scopeOfBioReplication = config$msstats$scopeOfBioReplication, scopeOfTechReplication = config$msstats$scopeOfTechReplication, interference = as.logical(config$msstats$interference), equalFeatureVar = as.logical(config$msstats$equalFeatureVar), missing.action = config$msstats$missing_action)$ComparisonResult
  results = groupComparison(data = mssquant, contrast.matrix = contrasts)
  write.table(results$ComparisonResult, file=config$files$output, eol="\n", sep="\t", quote=F, row.names=F, col.names=T)
  write.table(results$ModelQC, file=gsub(".txt","_ModelQC.txt",config$files$output), eol="\n", sep="\t", quote=F, row.names=F, col.names=T)
  cat(sprintf(">> WRITTEN\t%s\n",config$files$output))

  #(1) Minimal number of biological replicates per condition
  cat(">> CALCULATING SAMPLE SIZE FOR FUTURE EXPERIMENTS\n" )
  results.ss1 <- designSampleSize(data=results$fittedmodel,numSample=TRUE,desiredFC=c(1.25,2),FDR=0.05,power=0.95)
  results.ss2 <- designSampleSize(data=results$fittedmodel,numSample=TRUE,desiredFC=c(1.25,2),FDR=0.05,power=0.9)
  results.ss3 <- designSampleSize(data=results$fittedmodel,numSample=TRUE,desiredFC=c(1.25,2),FDR=0.05,power=0.8)
  results.ss = rbind( results.ss1, results.ss2, results.ss3)
  write.table(results.ss, file=gsub(".txt","_sampleSize.txt",config$files$output), eol="\n", sep="\t", quote=F, row.names=F, col.names=T)

  #(2) Power calculation
  cat(">> CALCULATING POWER OF EXPERIMENT\n" )
  results.power1 <- designSampleSize(data=results$fittedmodel,numSample=3, desiredFC=c(1.25,2),FDR=0.05,power=TRUE)
  results.power2 <- designSampleSize(data=results$fittedmodel,numSample=2, desiredFC=c(1.25,2),FDR=0.05,power=TRUE)
  results.power <- rbind(results.power1, results.power2)
  write.table(results.power, file=gsub(".txt","_experimentPower.txt",config$files$output), eol="\n", sep="\t", quote=F, row.names=F, col.names=T)

  return(results)
}


convertDataLongToMss = function(data_w, keys, config){
  cat(">> CONVERTING DATA TO MSSTATS FORMAT\n")
  data_l = meltMaxQToLong(data_w, na.rm = F)
  data_lk = mergeMaxQDataWithKeys(data_l, keys, by=c('RawFile','IsotopeLabelType'))
  dmss = dataToMSSFormat(data_lk)
  ## sanity check for zero's
  if(nrow(dmss[!is.na(dmss$Intensity) & dmss$Intensity == 0,]) > 0){
    dmss[!is.na(dmss$Intensity) & dmss$Intensity == 0,]$Intensity = NA
  }
  # Only write out if this was a data file. Otherwize, just pass the data back.
  if(!is.null(config$files$data)) write.table(dmss, file=gsub('.txt','-mss.txt',config$files$data), eol="\n", sep="\t", quote=F, row.names=F, col.names=T)
  return(dmss)
}


# annotate the files based on the Uniprot accession id's
Extras.annotate = function(results, output_file="MSstat_results.txt", uniprot_ac_col='Protein', group_sep=';', uniprot_dir = '~/github/kroganlab/source/db/', species='HUMAN'){
  cat(">> ANNOTATING\n")

  # remove unnamed proteins that are listed as ""
  if(length(which(results$Protein==""))>0)  results <- results[-which(results$Protein==""),]

  # read in all the annotation files from the uniprot_db directory
  species_split = unlist(strsplit(species, "-"))
  Uniprot = NULL
  for(org in species_split){
    cat(sprintf("\tLOADING %s\n",org))
    tmp = read.delim2(sprintf("%s/uniprot_protein_descriptions_%s.txt",uniprot_dir,org), stringsAsFactors=F, quote="")
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

# Annotate data, plot volcano plots, etc
writeExtras = function(results, config){

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
    heat_labels = prettyPrintHeatmapLabels(uniprot_acs=sign_hits$Protein,uniprot_ids=sign_hits$name, gene_names=sign_hits$gene_name)
    heat_data_w = plotHeat(mss_F = sign_hits, out_file =  gsub('.txt','-sign.pdf',config$files$output), names=heat_labels, cluster_cols=config$output_extras$heatmap_cluster_cols, display = config$output_extras$heatmap_display)
  }

  if(config$output_extras$volcano){
  	cat(">>   PLOTTING VOLCANO PLOT\n")
    file_name = gsub('.txt','-volcano.pdf',config$files$output)
    volcanoPlot(results_ann[grep(selected_labels,results_ann$Label),], lfc_upper, lfc_lower, FDR=config$output_extras$FDR, file_name=file_name)
  }
}

trim <- function (x) gsub("^\\s+|\\s+$", "", x)



#' @title Run MSStats pipeline designed for MaxQuant evidence files.
#' @description This function will take a MaxQuant file, keys, contrast, and config file, and handle all the file processing and conversions before running MSStats on the data. This is designed using the MSStats constraints/options outlined in the config.yaml file.
# @param input_file The filepath to the MaxQuant searched results file (txt tab delimited file). This can either be in either long or wide format.
# @param keys_file The filepath to the keys file used in the MSStats analysis.
# @param ref_proteome_file The filepath to the reference proteom (txt tab delimited file).
#' @param opt Pick which quantitative value to use: "spectral_counts" or "ms1"
#' @keywords MSStats MaxQuant
#' MQ.msstats()
MQ.msstats <- function(opt){
  opt$config = opt
  cat(">> MSSTATS PIPELINE\n")
  config = tryCatch(yaml.load_file(opt$config), error = function(e) {cat(opt$config);break} )

  # process MaxQuant data, link with keys, and convert for MSStats format
  if(config$data$enabled){
    cat(">> LOADING DATA\n")
    #data = fread(config$files$data, stringsAsFactors=F, integer64 = 'double')  # Found more bugs in fread. Not worth the compormise in data integrity just to save time reading in data
    data <- read.delim(config$files$data, stringsAsFactors=F, sep='\t')
    data <- data.table(data)
    setnames(data, colnames(data),gsub('\\s','.',colnames(data)))
    keys = read.delim(config$files$keys, stringsAsFactors=F, sep='\t')
    keys <- data.table(keys)

    ## the following lines were added to integrate the Label with the Filename when using multiple labels (e.g. H/L)
    ## currently we leave this in because the MSstats discinction on labeltype doesn't work
    ## see ISSUES https://github.com/everschueren/RMSQ/issues/1

    tryCatch(setnames(data, 'Raw.file', 'RawFile'), error=function(e) cat('Raw.file not found\n'))
    tryCatch(setnames(keys, 'Raw.file', 'RawFile'), error=function(e) cat('Raw.file not found\n'))

    cat('\tVERIFYING DATA AND KEYS\n')
    if(!'IsotopeLabelType' %in% colnames(data)) data[,IsotopeLabelType:='L']

    data = mergeMaxQDataWithKeys(data, keys, by = c('RawFile','IsotopeLabelType'))
    data$RawFile = paste(data$RawFile, data$IsotopeLabelType, sep='')
    keys$RawFile = paste(keys$RawFile, keys$IsotopeLabelType, sep='')
    keys$Run = paste(keys$IsotopeLabelType,keys$Run , sep='')
    data$IsotopeLabelType = 'L'
    keys$IsotopeLabelType = 'L'
    data[Intensity<1,]$Intensity=NA ## fix for weird converted values from fread

    ## end hacks for SILAC

    ## FILTERING : handles Protein Groups and Modifications
    if(config$filters$enabled) data_f = filterData(data, config) else data_f=data

    ## FORMATTING IN WIDE FORMAT FOR NORMALIZATION PURPOSES
    if(config$files$sequence_type == 'modified') castFun = castMaxQToWidePTM else castFun = castMaxQToWide
    data_w = castFun(data_f)

    ## test code for heatmap
    if(!is.null(config$files$sample_plots) && config$files$sample_plots){
      keys_in_data = keys[keys$RawFile %in% unique(data$RawFile),]
      sampleCorrelationHeatmap(data_w = data_w, keys = keys_in_data, config = config)
      samplePeptideBarplot(data_f, config)
    }

    ## AGGREGATION
    if(config$aggregation$enabled){
      res = aggregateData(data_w, keys, config)
      data_w=res$data_w_agg
      keys=res$keys_agg
    }

    ## NORMALIZATION -- we normalize with MSstats now
    #if(config$normalization$enabled) data_fn = normalizeData(data_w, config) else data_fn=data_w

  }

  ## MSSTATS
  if(config$msstats$enabled){

    if(is.null(config$msstats$msstats_input)){
      dmss = data.table(convertDataLongToMss(data_w, keys, config))

      ## make sure there are no doubles !!
      ## doubles could arise when protein groups are being kept and the same peptide is assigned to a unique Protein. Not sure how this is possible but it seems to be like this in maxquant output. A possible explanation is that these peptides have different retention times (needs to be looked into)
      dmss = data.frame(dmss[,j=list(ProteinName=paste(ProteinName,collapse=';'),Intensity=median(Intensity, na.rm=T)),by=c('PeptideSequence','ProductCharge','PrecursorCharge','FragmentIon','IsotopeLabelType','Run','BioReplicate','Condition')])

    }else{
      cat(sprintf("\tREADING PREPROCESSED\t%s\n", config$msstats$msstats_input))
      dmss = read.delim(config$msstats$msstats_input, stringsAsFactors=F, integer64 = 'double')
      dmss <- data.table(dmss)
    }

    cat(sprintf('>>LOADING MSSTATS %s VERSION\n', config$msstats$version))
    # check if there is nothing listed for the ms
    if( is.null(config$msstats$version) ){
      require('MSstats', character.only = T)
      cat('\tUsing Msstats version', installed.packages()['MSstats','Version'], '\n')
    }else if(config$msstats$version == 'MSstats.daily'){
      library('MSstats.daily', character.only = T)
    }else{
      stop( 'COULD NOT FIND MSSTATS VERSION. PLEASE MAKE SURE THERE IS A VERSION ON THIS MACHINE.')
    }

    # Read in contrasts file
    contrasts = read.delim(config$files$contrasts, stringsAsFactors=F)
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

#########################
## CONFIG LOADING #######

# spec = matrix(c(
#   'verbose', 'v', 2, "integer", "",
#   'help'   , 'h', 0, "logical", "available arguments (this screen)",
#   'config'  , 'c', 1, "character", "configuration file in YAML format"),
#   byrow=TRUE, ncol=5)
#
# opt = getopt(spec = spec, opt = commandArgs(TRUE), command = get_Rscript_filename(), usage = FALSE, debug = FALSE)
#
# # if help was asked for print a friendly message
# # and exit with a non-zero error code
# if ( !is.null(opt$help) ) {
#   cat(getopt(spec, usage=TRUE));
#   q(status=1);
# }

## TEST WORKS WITH LATEST CODE
# opt = c(opt, config='tests/LabelFree-ub/LabelFree-ub-test.yaml')

# if(!exists("DEBUG")){
#   cat(">> RUN MODE\n")
#   system.time(main(opt))
# }else{
#   cat(">> DEBUG MODE\n")
# }
# # rm(list = ls(all = TRUE))
