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
#' castMaxQToWide()
#' @export
castMaxQToWide <- function(d_long){
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
#' castMaxQToWidePTM()
#' @export
castMaxQToWidePTM <- function(d_long){
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
#' @title Filtering data
#' @description Apply the filtering options, i.e., remove protein groups and/or
#' contaminants, and/or, select posttranslational modification (if any)
#' @param data Evidence file (data.frame)
#' @param config Configuration object (opened yaml file)
#' @keywords filtering, remove, proteingroups, ptms
#' filterData()
#' @export

filterData <- function(data, config){
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

# ------------------------------------------------------------------------------
#' @title Remove contaminants and empty proteins
#' @description Remove contaminants and erronously identified 'reverse' 
#' sequences by MaxQuant
#' @param data the data.frame in long format
#' @return A new data.frame without REV__ and CON__ Proteins
#' @keywords cleanup, contaminants
#' filterMaxqData()
#' @export
filterMaxqData <- function(data){
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
#' mergeMaxQDataWithKeys()
#' @export
mergeMaxQDataWithKeys <- function(data, keys, by=c('RawFile')){
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
#' @title Convert the SILAC evidence file to MSstats format
#' 
#' @description Converting the evidence file from a SILAC search to a format 
#' compatible with MSstats. It basically modifies the Raw.files adding the 
#' Heavy and Light label
#' @param filename Text filepath to the evidence file
#' @param output Text filepath of the output name
#' @return A data.frame with SILAC data processed for MSstats
#' @keywords 
#' MQutil.SILACToLong()
#' @export
MQutil.SILACToLong = function(filename, output){
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
#' @title Read the Evidence File
#' @description Read in a MaxQuant searched Evidence file using data.table. 
#' This function properly classes each column and so fread doesn't have 
#' to guess.
#' @param evidence_file The filepath to the MaxQuant searched data (evidence) 
#' file (txt tab delimited file).
#' @return A data.frame with the evidence file with defining classes 
#' @keywords MaxQuant, evidence
#' read_evidence_file()
#' @export
read_evidence_file <- function(evidence_file){
  cat("Reading the evidence file...\n")
  # read in the first line to get the header names
  cols <- readLines(evidence_file, 1)
  cols <- data.frame( V1 = unlist(strsplit(cols, '\t')), stringsAsFactors = F)
  cols$idx <- 1:dim(cols)[1]
  
  # get data frame of pre-recorded column names and their respective classes
  col.classes <- as.data.frame( matrix(c("Sequence","character","Length","integer","Modifications","character","Modified sequence","character","Oxidation (M) Probabilities","character","Oxidation (M) Score Diffs","character","Acetyl (Protein N-term)","integer","Oxidation (M)","integer","Missed cleavages","integer","Proteins","character","Leading proteins","character","Leading razor protein","character","Gene names","character","Protein names","character","Type","character","Raw file","character","Experiment","character","MS/MS m/z","numeric","Charge","integer","m/z","numeric","Mass","numeric","Resolution","numeric","Uncalibrated - Calibrated m/z [ppm]","numeric","Uncalibrated - Calibrated m/z [Da]","numeric","Mass Error [ppm]","numeric","Mass Error [Da]","numeric","Uncalibrated Mass Error [ppm]","numeric","Uncalibrated Mass Error [Da]","numeric","Max intensity m/z 0","numeric","Retention time","numeric","Retention length","numeric","Calibrated retention time","numeric","Calibrated retention time start","numeric","Calibrated retention time finish","numeric","Retention time calibration","numeric","Match time difference","numeric","Match m/z difference","numeric","Match q-value","numeric","Match score","numeric","Number of data points","integer","Number of scans","integer","Number of isotopic peaks","integer","PIF","numeric","Fraction of total spectrum","numeric","Base peak fraction","numeric","PEP","numeric","MS/MS Count","integer","MS/MS Scan Number","integer","Score","numeric","Delta score","numeric","Combinatorics","integer","Intensity","numeric","Reverse","character","Potential contaminant","character","id","integer","Protein group IDs","character","Peptide ID","integer","Mod. peptide ID","integer","MS/MS IDs","character","Best MS/MS","integer","AIF MS/MS IDs","logical","Oxidation (M) site IDs","character", "Leading Proteins", "character", "Contaminant", "character"), ncol=2, byrow=T), stringsAsFactors = F)
  # merge the classes to the columns
  cols.matched = merge(cols, col.classes, by="V1", all.x=T)
  # re-order things to match the initial order
  cols.matched <- cols.matched[order(cols.matched$idx),]
  
  # Stop if there is an issue
  if(length(which(is.na(cols.matched$V2)))>0){
    stop(paste0("OH NO!! YOUR EVIDENCE FILE CONTAINS A COLUMN THAT I DON'T RECOGNIZE :( PLEASE TELL THE 'col.classes' IN THE read_evidence_file' FUNCTION AND ADD IN THIS NEW COLUMN(S) CALLED \n\t", paste(cols.matched$V1[which(is.na(cols.matched$V2))], collapse="\n\t"), "\n" ) )
  }
  
  # read in the evidence file with their classes
  x <- fread(evidence_file, integer64 = 'double', colClasses = cols.matched$V2)
  return(x)
}

# ------------------------------------------------------------------------------
#' @title Remove protein groups
#' @description Remove the group of proteins ids separated by separated by `;`
#' @param data Data.frame with a `Proteins` column.
#' @return A data.frame with the protein groups removed
#' @keywords maxquant, remove, proteingroups
#' removeMaxQProteinGroups()
#' @export
removeMaxQProteinGroups <- function(data){
  data_selected = data[grep(";",data$Proteins, invert=T),]
  return(data_selected)
}

# ------------------------------------------------------------------------------
#' @title Correlation heatmaps of all the individual features
#' @description Correlation heatmap using intensity values across all the 
#' conditions
#' @param data_w reshaped data.frame resulting from the `castMaxQToWidePTM` 
#' function
#' @param keys keys data.frame
#' @return A heatmap
#' @keywords heatmap, intensity, comparisons
#' sampleCorrelationHeatmap()
#' @export
sampleCorrelationHeatmap <- function (data_w, keys, config) {
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
#' samplePeptideBarplot()
#' @export
samplePeptideBarplot <- function(data_f, config){
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
#' significantHits()
#' @export
significantHits <- function(mss_results, labels='*', LFC=c(-2,2), FDR=0.05){
  ## get subset based on labels
  selected_results = mss_results[grep(labels,mss_results$Label), ]
  cat(sprintf('\tAVAILABLE LABELS FOR HEATMAP:\t%s\n',paste(unique(mss_results$Label), collapse=',')))
  cat(sprintf('\tSELECTED LABELS FOR HEATMAP:\t%s\n',paste(unique(selected_results$Label), collapse=',')))
  significant_proteins = selected_results[(!is.na(selected_results$log2FC) & selected_results$adj.pvalue <= FDR & (selected_results$log2FC >= LFC[2] | selected_results$log2FC <= LFC[1])) , 'Protein']
  significant_results = selected_results[selected_results$Protein %in% significant_proteins, ]
  return(significant_results)
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
#' @title Volcano plot (log2fc / pvalues)
#' @description It generates a scatter-plot used to quickly identify changes
#' @param mss_results_sel Selected MSstats results
#' @param lfc_upper log2fc upper threshold (positive value)
#' @param lfc_lower log2fc lower threshold (negative value)
#' @param FDR False Discovery Rate threshold
#' @param file_name Name for the output file
#' @param PDF Option to generate pdf format. Default: `T`
#' @param decimal_threshold Decimal threshold for the pvalue. 
#' Default: 16 (10^-16)
#' @keywords plot, volcano
#' volcanoPlot()
#' @export
volcanoPlot <- function(mss_results_sel, lfc_upper, lfc_lower, FDR, file_name='', PDF=T, decimal_threshold=16){
  
  # handle cases where log2FC is Inf. There are no pvalues or other information for these cases :(
  # Issues with extreme_val later if we have Inf/-Inf values.
  if( sum(is.infinite(mss_results_sel$log2FC)) > 0 ){
    idx <- is.infinite(mss_results_sel$log2FC)
    mss_results_sel$log2FC[ idx ] <- NA
  }
  
  min_x = -ceiling(max(abs(mss_results_sel$log2FC), na.rm=T))
  max_x = ceiling(max(abs(mss_results_sel$log2FC), na.rm=T))
  # Deal with special cases in the data where we have pvalues = Inf,NA,0
  if( sum(is.na(mss_results_sel$adj.pvalue))>0 ) mss_results_sel <- mss_results_sel[!is.na(mss_results_sel$adj.pvalue),]
  if(nrow(mss_results_sel[mss_results_sel$adj.pvalue == 0 | mss_results_sel$adj.pvalue == -Inf,]) > 0) mss_results_sel[!is.na(mss_results_sel$adj.pvalue) & (mss_results_sel$adj.pvalue == 0 | mss_results_sel$adj.pvalue == -Inf),]$adj.pvalue = 10^-decimal_threshold
  max_y = ceiling(-log10(min(mss_results_sel[mss_results_sel$adj.pvalue > 0,]$adj.pvalue, na.rm=T))) + 1
  
  l = length(unique(mss_results_sel$Label))
  w_base = 7
  h_base = 7
  
  if(l<=2){
    w=w_base*l 
  }else{
    w=w_base*2
  }
  h = h_base*ceiling(l/2)
  
  if(PDF) pdf(file_name, width=w, height=h)
  p = ggplot(mss_results_sel, aes(x=log2FC,y=-log10(adj.pvalue)))
  print(p + geom_point(colour='grey') + 
          geom_point(data = mss_results_sel[mss_results_sel$adj.pvalue <= FDR & mss_results_sel$log2FC>=lfc_upper,], aes(x=log2FC,y=-log10(adj.pvalue)), colour='red', size=2) +
          geom_point(data = mss_results_sel[mss_results_sel$adj.pvalue <= FDR & mss_results_sel$log2FC<=lfc_lower,], aes(x=log2FC,y=-log10(adj.pvalue)), colour='blue', size=2) +
          geom_vline(xintercept=c(lfc_lower,lfc_upper), lty='dashed') + 
          geom_hline(yintercept=-log10(FDR), lty='dashed') + 
          xlim(min_x,max_x) + 
          ylim(0,max_y) + 
          facet_wrap(facets = ~Label, ncol = 2, scales = 'fixed')) 
  if(PDF) dev.off()
}


# ------------------------------------------------------------------------------
#' @title Generate the contrast matrix required by MSstats from a txt file
#' @description It simplifies the process of creating the contrast file
#' @param contrast_file The text filepath of contrasts
#' @keywords 
#' writeContrast()
#' @export
writeContrast <- function(contrast_file) {
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