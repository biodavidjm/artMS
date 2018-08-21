# ------------------------------------------------------------------------------
#' @title Write extras
#' @description Extras after MSstats, as annotations, volcano plots, heatmaps
#' @param results MSstats results
#' @param config The configuration object (yaml)
#' @return Extras as selected in the yaml file, including:
#' - volcano plot (pdf)
#' - Adding annotations (gene symbol based on uniprot)
#' @keywords extras, annotations, volcano
#' artms_writeExtras()
#' @export
artms_writeExtras <- function(results, config){
  
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
  sign_hits <- artms_significantHits(results_ann,labels=selected_labels,LFC=c(lfc_lower,lfc_upper),FDR=config$output_extras$FDR)
  if( dim(sign_hits)[1] == 0 ) stop("NO SIGNIFICANT HITS DETECTED IN THIS EXPERIMENT. ABORTING PLOTS.\n")
  sign_labels <- unique(sign_hits$Label)
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
    artms_volcanoPlot(results_ann[grep(selected_labels,results_ann$Label),], lfc_upper, lfc_lower, FDR=config$output_extras$FDR, file_name=file_name)  
  }
}

# 
# ------------------------------------------------------------------------------
#' @title Annotate the files based on the Uniprot accession id
#' @description Annotate the files based on the Uniprot accession id
#' @param results MSstats results
#' @param output_file output file name
#' @param uniprot_ac_col Column with uniprot ids
#' @param group_sep Group separation character (if any. Eg: `;`)
#' @param uniprot_dir Directory with the Uniprot mappings
#' @param species Species (dash separated accordind to the uniprot file name)
#' @return Annotated data.frame
#' @keywords extras, annotations
#' Extras.annotate()
#' @export
Extras.annotate <- function(results, output_file, uniprot_ac_col='Protein', group_sep=';', uniprot_dir = '~/github/kroganlab/source/db/', species='HUMAN'){
  cat(">> ANNOTATING\n")
  
  # remove unnamed proteins that are listed as ""
  if(length(which(results$Protein==""))>0)  results <- results[-which(results$Protein==""),]
  
  # read in all the annotation files from the uniprot_db directory
  species_split = unlist(strsplit(species, "-"))
  Uniprot = NULL
  for(org in species_split){
    cat(sprintf("\tLOADING %s\n",org))
    tmp <- read.delim(sprintf("%s/uniprot_protein_descriptions_%s.txt",uniprot_dir,org), stringsAsFactors=F, quote="")    
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
