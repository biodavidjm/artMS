# ==============================================================================
## Small MSstats-related Functions

# ------------------------------------------------------------------------------
#' @title Long to Wide format using the `Sequence` column of the evidence file
#' 
#' @description Facilitates applying the dcast function, i.e., takes long-format
#' data and casts it into wide-format data.
#' @param d_long the data.frame in long format
#' @return Evidence file reshaped by rawfile and IsotopeLabelType
#' @keywords data.frame, dcast
#' artms_castMaxQToWide()
#' @export
artms_castMaxQToWide <- function(d_long){
  data_w = data.table::dcast( Proteins + Sequence + Charge ~ RawFile + IsotopeLabelType, data=d_long, value.var='Intensity', fun.aggregate=sum, fill = NA)
  return(data_w)
}

# ------------------------------------------------------------------------------
#' @title Long to Wide format selecting the `Modified.sequence` column of the 
#' evidence file
#' @description Facilitates applying the dcast function, i.e., takes long-format 
#' data and casts it into wide-format data.
#' @param d_long the data.frame in long format
#' @return Evidence file reshaped by rawfile and IsotopeLabelType
#' @keywords data.frame, dcast
#' artms_castMaxQToWidePTM()
#' @export
artms_castMaxQToWidePTM <- function(d_long){
  data_w = data.table::dcast( Proteins + Modified.sequence + Charge ~ RawFile + IsotopeLabelType, data=d_long, value.var='Intensity', fun.aggregate=sum, fill=NA)
  setnames(data_w,2,'Sequence')
  return(data_w)
}

# ------------------------------------------------------------------------------
#' @title Change a specific column name in a given data.frame
#' @description Making easier j
#' @param dataset the data.frame with the column name you want to change
#' @param oldname the old column name
#' @param neename the new name for that column
#' @return A new data.frame with the new specified column name
#' @keywords rename, data.frame, columns
#' changeColumnName()
#' @export
changeColumnName <- function(dataset, oldname, newname){
  if( !(oldname %in% colnames(dataset)) ){
    stop("The Column name provided <",oldname,"> was not found in the data.table provided")
  }
  names(dataset)[grep(paste0('^',oldname,'$'), names(dataset))] <- newname
  return(dataset)
}

# ------------------------------------------------------------------------------
#' @title Protein abundance dot plots
#' 
#' @description Protein abundance dot plots for each unique uniprot id
#' @param input_file The `-normalized.txt` output file from MSstats
#' @param output_file Wished output file name (add the `.pdf` extension)
#' @return A pdf file with each individual protein abundance plot for each
#' conditions
#' @keywords abundance, dotplots, plot
#' artms_dataPlots()
#' @export
artms_dataPlots <- function(input_file, output_file){
  
  data_mss = fread(input_file, integer64 = 'double')
  unique_subjects = unique(data_mss$PROTEIN)
  condition_length = length(unique(data_mss$GROUP_ORIGINAL))
  min_abu = min(data_mss$ABUNDANCE, na.rm = T)
  max_abu = max(data_mss$ABUNDANCE, na.rm=T)
  
  pdf(output_file, width = condition_length*1.5, height = 3)
  
  cat('PRINTING CONDITION PLOTS\n')
  for(subject in unique_subjects){
    subject_data = data_mss[PROTEIN==subject,]
    cat(sprintf('\t%s\n',subject))
    p = ggplot(data = subject_data, aes(x=SUBJECT_ORIGINAL,y=ABUNDANCE, colour=FEATURE))
    p = p + geom_point(size=2) + 
      facet_wrap(facets = ~ GROUP_ORIGINAL, drop = T, scales = 'free_x', ncol = condition_length) + 
      ylim(min_abu,max_abu) +
      theme(axis.text.x=element_text(angle=-90,hjust=1)) +
      guides(colour=FALSE) +
      xlab(NULL) +
      ggtitle(subject)
    print(p)
  }
  dev.off()
}

# ------------------------------------------------------------------------------
#' @title Filtering data
#' @description Apply the filtering options, i.e., remove protein groups and/or
#' contaminants, and/or, select posttranslational modification (if any)
#' @param data Evidence file (data.frame)
#' @param config Configuration object (opened yaml file)
#' @keywords filtering, remove, proteingroups, ptms
#' artms_filterData()
#' @export
artms_filterData <- function(data, config){
  cat(">> FILTERING\n")
  if(config$filters$protein_groups == 'remove'){
    cat("\tPROTEIN GROUPS\tREMOVE\n")
    data_f = removeMaxQProteinGroups(data)  
  }else if(config$filters$protein_groups == 'keep'){
    cat("\tPROTEIN GROUPS\tIGNORE\n")
    data_f = data
  }else{
    stop("\n\nFILTERING OPTION FOR protein_groups NOT UNDERSTOOD (OPTIONS AVAILABLE: keep OR remove\n\n")
  }
  
  if(config$filters$contaminants){
    cat("\tCONTAMINANTS\tREMOVE\n")
    data_f = artms_filterMaxqData(data_f)  
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

# ------------------------------------------------------------------------------
#' @title Remove contaminants and empty proteins
#' @description Remove contaminants and erronously identified 'reverse' 
#' sequences by MaxQuant
#' @param data the data.frame in long format
#' @return A new data.frame without REV__ and CON__ Protein ids
#' @keywords cleanup, contaminants
#' artms_filterMaxqData()
#' @export
artms_filterMaxqData <- function(data){
  # Remove contaminants and reversed sequences (labeled by MaxQuant)
  data_selected = data[grep("CON__|REV__",data$Proteins, invert=T),]
  # Remove empty proteins names
  blank.idx <- which(data_selected$Proteins == "")
  if(length(blank.idx)>0)  data_selected = data_selected[-blank.idx,]
  return(data_selected)
}

# ------------------------------------------------------------------------------
#' @title Merge evidence and keys files
#' @description Merge the evidence and keys files on the given columns
#' @param data The evidence in data.frame
#' @param keys The keys in data.frame
#' @param by Vector specifying the columns use to merge the evidence and keys.
#' @return A new data.frame with the evidence and keys merged
#' Default: `RawFile`
#' @keywords merge, evidence, keys
#' artms_mergeMaxQDataWithKeys()
#' @export
artms_mergeMaxQDataWithKeys <- function(data, keys, by=c('RawFile')){
  # Check if the number of RawFiles is the same.
  unique_data <- unique(data$RawFile)
  unique_keys <- unique(keys$RawFile)
  
  if (length(unique_keys) != length(unique_data)){
    keys_not_found = setdiff(unique_keys, unique_data)
    data_not_found = setdiff(unique_data, unique_keys)
    cat(sprintf("\tkeys found: %s \n\t keys not in data file:\n%s\n", length(unique_keys)-length(keys_not_found), paste(keys_not_found,collapse='\t')))
    cat(sprintf("\tdata found: %s \n\t data not in keys file:\n%s\n", length(unique_data)-length(data_not_found), paste(data_not_found, collapse='\t')))
  }else{
    cat("\t>-----+ Check point: the number of RawFiles in both keys and evidences file is identical\n\n")
  }
  
  ## select only required attributes from MQ format
  data = merge(data, keys, by=by)
  return(data)
}

# ------------------------------------------------------------------------------
#' @title Merge evidence and keys by file name
#' 
#' @description Merge evidence and keys by file name
#' @param evidence_file The Evidence file name
#' @param keys_file The keys file name
#' @return A data.frame with both evidence and keys files merged by raw.files
#' @keywords merge, evidence, keys
#' artms_mergeEvidenceKeysByFiles()
#' @export
artms_mergeEvidenceKeysByFiles <- function(evidence_file, keys_file) {
  
  data <- read.delim(evidence_file, sep='\t', quote = "", header = T, stringsAsFactors = F)
  keys <- read.delim(keys_file, sep='\t', quote = "", header = T, stringsAsFactors = F)
  
  if( !('RawFile' %in% colnames(data)) ) {
    tryCatch(changeColumnName(data, 'Raw.file', 'RawFile'), error=function(e) stop('\nRaw.file not found in EVIDENCE FILE\n'))
  }
  
  if( !('RawFile' %in% colnames(keys)) ) {
    tryCatch(changeColumnName(keys, 'Raw.file', 'RawFile'), error=function(e) stop('\nRaw.file not found in EVIDENCE FILE\n'))
  }
  
  # Check that the keys file is correct
  if(any(!c('RawFile','IsotopeLabelType','Condition','BioReplicate','Run') %in% colnames(keys))){ #,'SAINT','BioReplicaSaint'
    cat('\n-!-!-!-!- ERROR: COLUMN NAMES IN KEYS NOT CONFORM TO SCHEMA. One of these is lost\n
        \tRawFile\n\tIsotopeLabelType\n\tCondition\n\tBioReplicate\n\tRun\n') # \tSAINT\n\tBioReplicaSaint\n\n
    stop('Please, try again once revised\n\n')
  }
  
  # MERGING THE DATA
  # Checking that the keys make sense
  unique_data <- unique(data$RawFile)
  unique_keys <- unique(keys$RawFile)
  
  # Rawfiles on Keys not found on the data
  keys_not_found <- setdiff(unique_keys, unique_data)
  # Rawfiles on Data not found on the keys
  data_not_found <- setdiff(unique_data, unique_keys)
  
  if ( (length(keys_not_found) != 0) & ( length(data_not_found) != 0) ) {
    cat(sprintf("keys found: %s \t keys not in data file:\n%s\n", length(unique_keys)-length(keys_not_found), paste(keys_not_found,collapse='\t')))
    cat(sprintf("data found: %s \t data not in keys file:\n%s\n", length(unique_data)-length(data_not_found), paste(data_not_found, collapse='\t')))
    stop('\nThis script is sorry, but it needs to stop this because something is going on between your keys and evidence files so you better check\n')
  }
  
  ## select only required attributes from MQ format
  datamerged <- merge(data, keys, by='RawFile')
  
  return(datamerged)
}


# ------------------------------------------------------------------------------
#' @title Convert the SILAC evidence file to MSstats format
#' 
#' @description Converting the evidence file from a SILAC search to a format 
#' compatible with MSstats. It basically modifies the Raw.files adding the 
#' Heavy and Light label
#' @param filename Text filepath to the evidence file
#' @param output Text filepath of the output name
#' @return A data.frame with SILAC data processed for MSstats
#' @keywords 
#' artms_SILACtoLong()
#' @export
artms_SILACtoLong <- function(filename, output){
  file = Sys.glob(filename)
  cat(sprintf('>> PROCESSING SILAC EVIDENCE FILE\n'))
  tmp = fread(file, integer64 = 'double')
  
  # reshape the data and split the heavy and light data
  tmp_long = reshape2::melt(tmp, measure.vars = c('Intensity L','Intensity H'))
  tmp_long[,Intensity:=NULL]
  setnames(tmp_long,'value','Intensity')
  setnames(tmp_long,'variable','IsotopeLabelType')
  setnames(tmp_long,'Raw file','Raw.file')
  levels(tmp_long$IsotopeLabelType) = c('L','H')
  tmp_long[!is.na(tmp_long$Intensity) && tmp_long$Intensity<1,]$Intensity=NA
  write.table(tmp_long, file=output, sep='\t', quote=F, row.names=F, col.names=T)
  cat("----- + File ",output, " has been created\n")
  return(tmp_long)
}

# ------------------------------------------------------------------------------
#' @title Pretty Labels for Heatmaps
#' @description Generates pretty labels for the heatmaps.
#' @param uniprot_acs Uniprot accession id
#' @param uniprot_ids Uniprot entry id
#' @param gene_names Gene symbol
#' @return Pretty labels for a heatmap
#' @keywords plots, pretty
#' prettyPrintHeatmapLabels()
#' @export
prettyPrintHeatmapLabels <- function(uniprot_acs, uniprot_ids, gene_names){
  result = paste(uniprot_acs,uniprot_ids,gene_names,sep=' ')
  return(result)
}

# ------------------------------------------------------------------------------
#' @title Reshape the MSstats results file from long to wide format
#' 
#' @description Converts the normal MSStats output file into "wide" format 
#' where each row represents a protein's results, and each column represents 
#' the comparison made by MSStats. The fold change and p-value of each 
#' comparison will be it's own column.
#' @param input_file Input file name and location (MSstats `results.txt` file)
#' @param output_file Output file name and location
#' @return Reshaped file with unique protein ids and as many columns log2fc
#' and adj.pvalues as comparisons available
#' for as many 
#' @keywords
#' artms_resultsWide()
#' @export
artms_resultsWide <- function(input_file, output_file){
  input = fread(input_file, integer64 = 'double')
  input_l = melt(data = input[,c('Protein', 'Label','log2FC','adj.pvalue'), with=F], id.vars = c('Protein', 'Label'))
  
  ## then cast to get combinations of LFCV/PVAl and Label as columns
  input_w = dcast.data.table( Protein ~ Label+variable, data=input_l, value.var=c('value'))
  write.table(input_w, file=output_file, eol='\n', sep='\t', quote=F, row.names=F, col.names=T)
}

# ------------------------------------------------------------------------------
#' @title Correlation heatmaps of all the individual features
#' @description Correlation heatmap using intensity values across all the 
#' conditions
#' @param data_w reshaped data.frame resulting from the `artms_castMaxQToWidePTM` 
#' function
#' @param keys keys data.frame
#' @param config Configuration object (yaml loaded)
#' @return A heatmap
#' @keywords heatmap, intensity, comparisons
#' artms_sampleCorrelationHeatmap()
#' @export
artms_sampleCorrelationHeatmap <- function (data_w, keys, config) {
  mat = log2(data_w[,4:ncol(data_w),with=F])
  mat[is.na(mat)]=0
  mat_cor = cor(mat, method = 'pearson', use = 'everything')
  ordered_keys = keys[with(keys, order(RawFile)),] ## we want to make informarive row names so order by RawFile because that's how data_w is ordered
  mat_names = paste(ordered_keys$Condition, ordered_keys$BioReplicate, ordered_keys$Run)
  colnames(mat_cor) = mat_names
  rownames(mat_cor) = mat_names
  colors_pos = colorRampPalette(RColorBrewer::brewer.pal("Blues",n=5))(10)
  colors_neg = rev(colorRampPalette(RColorBrewer::brewer.pal("Reds",n=5))(10))
  colors_tot = c(colors_neg, colors_pos)
  pheatmap(mat = mat_cor, cellwidth = 10, cellheight = 10, scale = 'none', filename = gsub('.txt','-heatmap.pdf',config$files$output), color = colors_tot, breaks = seq(from=-1,to = 1, by=.1), fontfamily="mono")
  
}

# ------------------------------------------------------------------------------
#' @title Barplot of peptide counts per biological replicate
#' @description Total number of unique peptide identified per biological 
#' replicate
#' @param data_f Evidence file (same structure as the original)
#' @param config Configuration object
#' @return Barplot of peptide counts
#' @keywords barplot, counts, peptides
#' artms_samplePeptideBarplot()
#' @export
artms_samplePeptideBarplot <- function(data_f, config){
  # set up data into ggplot compatible format
  data_f = data.table(data_f, labels=paste(data_f$RawFile, data_f$Condition, data_f$BioReplicate))
  data_f = data_f[with(data_f, order(labels,decreasing = T)),]
  
  # plot the peptide counts for all the samples TOGETHER
  p = ggplot(data = data_f, aes(x=labels))
  p = p + geom_bar() + theme(axis.text.x = element_text(angle = 90, hjust = 1, family = 'mono')) + ggtitle('Unique peptides per run\n after filtering') + coord_flip()
  ggsave(filename = gsub('.txt','-peptidecounts.pdf',config$files$output), plot=p, width = 8, height = 10)
  
  w = 10
  h = ceiling( (7/5+2) * ceiling(length(unique(data_f$Condition))/5) )
  # plot the peptide counts for all the samples PER BAIT
  p = ggplot(data = data_f, aes(x=as.factor(BioReplicate)))
  p = p + geom_bar() + theme(axis.text.x = element_text(angle = 90, hjust = 1, family = 'mono')) + ggtitle('Unique peptides per run\n after filtering') + facet_wrap(~Condition, scales='free', ncol=5)  + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(filename = gsub('.txt','-peptidecounts-perBait.pdf',config$files$output), plot=p, width = w, height = h)
  
}

# ------------------------------------------------------------------------------
#' @title Select significant hits
#' @description Filtered data.frame with significant values (log2fc > 2 | 
#' log2fc < -2; adj.pvalue < 0.05)
#' @param mss_results data.frame of MSstats results
#' @param labels vector of selected labels. Default: all (`*`)
#' @param LFC log2fc vector, with the negative and positive threshold. Default: 
#' c(-2, 2)
#' @param FDR false discovery rate (adj.pvalue) threshold. Default: 0.05
#' @return A data.frame only with significant hits
#' @keywords significant, selections
#' artms_significantHits()
#' @export
artms_significantHits <- function(mss_results, labels='*', LFC=c(-2,2), FDR=0.05){
  ## get subset based on labels
  selected_results = mss_results[grep(labels,mss_results$Label), ]
  cat(sprintf('>> AVAILABLE LABELS FOR HEATMAP: %s\n',paste(unique(mss_results$Label), collapse=', ')))
  cat(sprintf('>> SELECTED LABELS FOR HEATMAP: %s\n',paste(unique(selected_results$Label), collapse=', ')))
  significant_proteins = selected_results[(!is.na(selected_results$log2FC) & selected_results$adj.pvalue <= FDR & (selected_results$log2FC >= LFC[2] | selected_results$log2FC <= LFC[1])) , 'Protein']
  significant_results = selected_results[selected_results$Protein %in% significant_proteins, ]
  return(significant_results)
}


# ------------------------------------------------------------------------------
#' @title Outputs the spectral counts from the MaxQuant evidence file.
#' 
#' @description Outputs the spectral counts from the MaxQuant evidence file.
#' @param input_file Maxquant evidence file
#' @param keys_file Keys file with the experimental design
#' @param output_file Output file name (add `.txt` extenstion)
#' @return A txt file with biological replicates, protein id, and spectral 
#' count columns
#' @keywords spectral_counts, evidence
#' artms_spectralCounts()
#' @export
artms_spectralCounts <- function(input_file, keys_file, output_file){
  data = fread(input_file, integer64 = 'double')
  keys = fread(keys_file, integer64 = 'double')
  
  tryCatch(setnames(data, 'Raw file', 'RawFile'), error=function(e) cat('Raw file not found. Try searching Raw.file instead\n'))
  tryCatch(setnames(keys, 'Raw.file', 'RawFile'), error=function(e) cat('Raw file not found. Try searching Raw file instead\n'))
  
  cat('\tVERIFYING DATA AND KEYS\n')
  if(!'IsotopeLabelType' %in% colnames(data)) data[,IsotopeLabelType:='L']
  data <- artms_mergeMaxQDataWithKeys(data, keys, by = c('RawFile','IsotopeLabelType'))
  data_sel <- data[,c('Proteins', 'Condition', 'BioReplicate', 'Run', 'MS/MS Count'), with=F]
  setnames(data_sel, 'MS/MS Count', 'spectral_counts')
  data_sel = aggregate( spectral_counts ~ Proteins+Condition+BioReplicate+Run, data=data_sel, FUN = sum)
  data_sel = data.frame(data_sel, bait_name=paste(data_sel$Condition, data_sel$BioReplicate, data_sel$Run, sep='_'))
  write.table(data_sel[,c('bait_name','Proteins','spectral_counts')], file=output_file, eol='\n', sep='\t', quote=F, row.names=F, col.names=T)
}

# ------------------------------------------------------------------------------
#' @title Remove white spaces
#' @description Remove white spaces
#' @param x A string
#' @keywords remove, whitespace
#' trim()
#' @export
trim <- function (x){
  gsub("^\\s+|\\s+$", "", x)
}


# ------------------------------------------------------------------------------
#' @title Generate the contrast matrix required by MSstats from a txt file
#' @description It simplifies the process of creating the contrast file
#' @param contrast_file The text filepath of contrasts
#' @keywords 
#' artms_writeContrast()
#' @export
artms_writeContrast <- function(contrast_file) {
  input_contrasts <- readLines(contrast_file, warn=F)
  
  # check if contrast_file is old-format (i.e the contrast_file is a matrix)
  headers <- unlist(strsplit(input_contrasts[1], split = "\t"))
  if (length(headers)>1) {
    newinput_contrasts <- c()
    for (i in 2:length(input_contrasts)){
      newinput_contrasts <- c(newinput_contrasts, unlist(strsplit(input_contrasts[i], split = "\t"))[1])
    }
    input_contrasts <- newinput_contrasts
  }
  
  # validate the input
  input_contrasts <- trimws(input_contrasts)
  valid <- TRUE
  accepted_chars <- c(LETTERS, letters, 0:9, '-','_')
  for (x in input_contrasts) {
    characs <- unlist(strsplit(x, split=''))
    not_allowed_count <- length(which(!(characs %in% accepted_chars)))
    if (not_allowed_count > 0) {
      valid <- FALSE
      cat(paste(x, "is not a valid input"))
      break
    }
    
    dash_count <- length(which(characs == '-'))
    if (dash_count != 1) {
      valid <- FALSE
      cat(paste(x, "needs to contain exactly 1 '-'"))
      break    
    }
  }
  
  if (valid) {
    mat <- t(as.data.frame(strsplit(input_contrasts, split = '-')))
    rownames(mat) <- NULL
    conds <- sort(unique(c(mat[,1], mat[,2])))
    contrast_matrix <- matrix(0, nrow = nrow(mat), ncol = length(conds))
    colnames(contrast_matrix) <- conds
    rownames(contrast_matrix) <- input_contrasts
    
    for (i in 1:nrow(mat)) {
      cond1 <- mat[i,1]
      cond2 <- mat[i,2]
      contrast_matrix[i, cond1] <- 1
      contrast_matrix[i, cond2] <- -1
    }
    return (contrast_matrix)
  }
  return (NA)
}
