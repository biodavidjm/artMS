# ------------------------------------------------------------------------------
#' @title Quality Control analysis of the MaxQuant evidence file
#' 
#' @description Quality Control analysis of the MaxQuant evidence file
#' @param evidence_file The evidence file
#' @param keys_file The keys file with the experimental design
#' @param prot_exp Proteomics experiment. 4 options available:
#' - `apms`: affinity purification mass spectrometry
#' - `ab`: protein abundance
#' - `ph`: protein phosphorylation
#' - `ub`: protein ubiquitination (aka ubiquitylation)
#' @param output_dir output directory to save results in
#' @param fractions Is a fractionated experiment?
#' - 1 yes
#' - 0 no (default)
#' @return Quality control files and plots
#' @keywords QC, quality, control, evidence
#' artms_evidenceQC()
#' @export
artms_evidenceQC.R <- function(evidence_file, keys_file, prot_exp, output_dir, fractions){
  cat("\n\n\t~~> Checking the Evidence file\n")
  cat("\t(it should take some time due to the usual large size of evidence files)\n")
  
  # EVIDENCE:
  evidencekeys <- artms_mergeEvidenceKeysByFiles(evidence_file, keys_file)
  
  if(fractions){
    # Check that the keys file is correct
    keys <- read.delim(keys_file, sep='\t', quote = "", header = T, stringsAsFactors = F)
    if(any(!c('RawFile','IsotopeLabelType','Condition','BioReplicate','Run', 'FractionKey') %in% colnames(keys))){
      cat('\nERROR: COLUMN NAMES IN KEYS NOT CONFORM TO SCHEMA. One of these is lost\n
          \tRawFile\n\tIsotopeLabelType\n\tCondition\n\tBioReplicate\n\tRun\n\n\tFractionKey\n\n\tSAINT\n\n')
      stop('Please, try again once revised\n\n')
    }
  }
  
  ekselecta <- aggregate(Intensity~Proteins+Condition+BioReplicate+Run, data=evidencekeys, FUN = sum)
  ekselectaBioreplica <- aggregate(Intensity~Proteins+Condition+BioReplicate, data=ekselecta, FUN = sum) 
  
  # create output directory if it doesn't exist
  if(!dir.exists(output_dir)){
    dir.create(output_dir, recursive = T)
  }
  
  # Checking the overall distribution of intensities before anything else
  # Based on Intensity
  ekselectaBioreplica$Abundance <- log2(ekselectaBioreplica$Intensity)
  j <- ggplot(ekselectaBioreplica, aes(BioReplicate,Intensity))
  j <- j + geom_jitter(width = 0.3, size = 0.5)
  j <- j + theme_minimal()
  j <- j + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
  # Based on Abundance
  k <- ggplot(ekselectaBioreplica, aes(BioReplicate,Abundance))
  k <- k + geom_jitter(width = 0.3, size = 0.5)
  k <- k + theme_minimal()
  k <- k + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
  
  cat("\t~~> Printing out the intensity distribution plots\n")
  intDistribution <- gsub("evidence.txt", "IntensityDistributions.pdf", evidence_file)
  intDistribution <- paste0("plot.",prot_exp,".",intDistribution)
  intDistribution <- paste0(output_dir,"/",intDistribution)
  
  pdf(intDistribution)
    plot(j)
    plot(k)
  garbage <- dev.off()
  
  # Feature generation: Combine Sequence and Charge.
  evidencekeys$Feature <- paste0(evidencekeys$Modified.sequence,"_",evidencekeys$Charge)
  
  # ============================================================================
  # GENERAL QUALITY CONTROL: CHECK PROPORTION OF CONTAMINANTS
  # sum of total intensities, label them as contaminants and non-contaminants
  # and plot intensities for each group
  
  # Careful with old versions of MaxQuant:
  if( any(grep("Leading.Proteins", names(evidencekeys))) ){
    evidencekeys <- changeTheName(evidencekeys,"Leading.Proteins","Leading.proteins")
  }
  
  # Combine all the fractions if this is a fractioning experiment by summing 
  # them up
  if (fractions){
    # Sum up all the fractions first
    evidencekeys <- aggregate(Intensity~Feature+Proteins+Leading.proteins+Condition+BioReplicate+Run, data=evidencekeys, FUN = sum)
  }
  
  evigeneral <- evidencekeys[c('Feature', 'Proteins', 'Leading.proteins', 'Intensity', 'Condition', 'BioReplicate','Run')]
  
  # Now let's classify the proteins as contaminants and no contaminants
  evigeneral$Contaminant <- ifelse(grepl("CON_|REV_",evigeneral$Leading.proteins), ifelse(grepl("REV_",evigeneral$Leading.proteins),"REV","CON"), "PROT")
  
  # PRINT OUT STATS about contaminants
  numberContaminants <- aggregate(data = evigeneral, Leading.proteins~Contaminant, FUN = length)
  print(numberContaminants)
  
  # AGGREGATE ALL THE INTENSITIES PER PROTEIN, summing everything up
  evigeneral$Intensity <- as.numeric(evigeneral$Intensity)
  
  evisummary <- aggregate(Intensity~Condition+BioReplicate+Contaminant+Run, data=evigeneral, FUN = sum)
  
  # CLEANING THE EVIDENCE OF CONTAMINANTS
  evidencekeysclean <- artms_filterMaxqData(evidencekeys)
  cat("\t~~> Data loaded\n")
  
  cat("\t~~> Printing out the reproducibility plots\n")
  seqReproName <- gsub("evidence.txt", "basicReproducibility.pdf", evidence_file)
  seqReproName <- paste0("plot.",prot_exp,".",seqReproName)
  seqReproName <- paste0(output_dir,"/",seqReproName)
  
  if (prot_exp == "ub") {
    evidencekeysclean <- evidencekeysclean[c('Feature', 'Modified.sequence', 'Proteins', 'Intensity', 'Condition', 'BioReplicate', 'Run')]
    evidencekeysclean <- evidencekeysclean[grep("(gl)", evidencekeysclean$Modified.sequence),]  
  }else if (prot_exp == "ph" | prot_exp == "ptmph") {
    evidencekeysclean <- evidencekeysclean[c('Feature', 'Modified.sequence', 'Proteins', 'Intensity', 'Condition', 'BioReplicate', 'Run')]
    evidencekeysclean <- evidencekeysclean[grep("(ph)", evidencekeysclean$Modified.sequence),]
  }else{
    cat("Good to go!\n")
  }
  
  pdf(seqReproName)
    artms_plotReproducibilityEvidence(evidencekeysclean)
  garbage <- dev.off()
  
  
  # Create matrix of reproducibility TECHNICAL REPLICAS
  data2matrix <- evidencekeysclean
  # Make sure that the Intensity is numeric
  data2matrix$Intensity <- as.numeric(data2matrix$Intensity)
  
  # Check the number of TECHNICAL REPLICAS by checking the first technical replica
  technicalReplicas <- unique(data2matrix$Run[which(data2matrix$BioReplicate == data2matrix$BioReplicate[1])])
  palette.breaks <- seq(1, 3, 0.1)
  color.palette  <- colorRampPalette(c("white","steelblue"))(length(palette.breaks))
  
  
  if(length(technicalReplicas) > 1){
    # First aggregate at the protein level by summing up everything
    biorepliaggregated <- aggregate(Intensity~Feature+Proteins+Condition+BioReplicate+Run, data = data2matrix, FUN = sum)
    biorepliaggregated$Intensity <- log2(biorepliaggregated$Intensity)
    evidencekeyscleanDCASTbioreplicas <- dcast(data=biorepliaggregated, Proteins+Feature~BioReplicate+Run, value.var = "Intensity")
    precordfBioreplicas <- evidencekeyscleanDCASTbioreplicas[,3:dim(evidencekeyscleanDCASTbioreplicas)[2]]
    Mtechnicalrep <- cor(precordfBioreplicas, use = "pairwise.complete.obs")
    
    theTechCorDis <- artms_plotCorrelationDistribution(Mtechnicalrep)
    
    # And now for clustering
    cat("\t~~> Correlation matrix including Technical replicas\n")
    matrixCorrelationBioreplicas <- gsub("evidence.txt", "correlationMatrixTR.pdf", evidence_file)
    matrixCorrelationBioreplicas <- paste0("plot.",prot_exp,".",matrixCorrelationBioreplicas)
    matrixCorrelationBioreplicas <- paste0(output_dir,"/",matrixCorrelationBioreplicas)
    
    pdf(matrixCorrelationBioreplicas, width = 20, height = 20) #
    corrplot(Mtechnicalrep, 
             method="square", 
             addCoef.col = "white",
             number.cex=1, 
             tl.cex = 1.5,
             tl.col = "black", 
             title = "Matrix Correlation based on peptide Intensities")
    pheatmap(Mtechnicalrep, 
             cluster_rows = T,
             cluster_cols = T,
             cellheight = 10,
             cellwidth=25, 
             main = "Clustering Technical Replicates",
             fontsize=6, 
             fontsize_row=8, 
             fontsize_col=12, 
             border_color='black',
             fontfamily="Helvetica",
             treeheight_row = F, 
             treeheight_col = F,
             color = color.palette
    )
    print(theTechCorDis)
    garbage <- dev.off()
  }else{
    cat("---NO TECHNICAL REPLICATES DETECTED\n")
  }

  # biological replicates  
  biorepliaggregated <- NULL

  # First deal with the technical replica
  if(length(technicalReplicas) > 1){
    # Aggregate at the protein level by summing up everything
    biorepliaggregated <- aggregate(Intensity~Feature+Proteins+Condition+BioReplicate+Run, data = data2matrix, FUN = sum)
    biorepliaggregated <- aggregate(Intensity~Feature+Proteins+Condition+BioReplicate, data = biorepliaggregated, FUN = mean)
  }else{
    biorepliaggregated <- aggregate(Intensity~Feature+Proteins+Condition+BioReplicate, data = data2matrix, FUN = sum)
  }
  
  biorepliaggregated$Intensity <- log2(biorepliaggregated$Intensity)
  evidencekeyscleanDCASTbioreplicas <- dcast(data=biorepliaggregated, Proteins+Feature~BioReplicate, value.var = "Intensity")
  precordfBioreplicas <- evidencekeyscleanDCASTbioreplicas[,3:dim(evidencekeyscleanDCASTbioreplicas)[2]]
  Mbioreplicas <- cor(precordfBioreplicas, use = "pairwise.complete.obs")
  
  theBiorCorDis <- artms_plotCorrelationDistribution(Mbioreplicas)
  
  cat("\t~~> BioReplicate correlation matrix (based on Proteins)\n")
  matrixCorrelationBioreplicas <- gsub("evidence.txt", "correlationMatrixBR.pdf", evidence_file)
  matrixCorrelationBioreplicas <- paste0("plot.",prot_exp,".",matrixCorrelationBioreplicas)
  matrixCorrelationBioreplicas <- paste0(output_dir,"/",matrixCorrelationBioreplicas)
  pdf(matrixCorrelationBioreplicas, width = 20, height = 20) #, width = 20, height = 20
  corrplot(Mbioreplicas, 
           method="square", 
           addCoef.col = "white",
           number.cex=0.9,
           tl.cex = 1.5, tl.col = "black", title = "Matrix Correlation based on Protein Intensities")
  pheatmap(Mbioreplicas, 
           cluster_rows = T,
           cluster_cols = T,
           cellheight = 10, 
           cellwidth=25, 
           main = "Clustering Biological Replicates",
           fontsize=6, 
           fontsize_row=8, 
           fontsize_col=12, 
           border_color='black',
           fontfamily="Helvetica",
           treeheight_row = F, 
           treeheight_col = F,
           color = color.palette
  )
  print(theBiorCorDis)
  garbage <- dev.off()
  
  # Create matrix of reproducibility CONDITIONS
  biorepliaggregated <- NULL
  
  # biological replicas
  # First deal with the technical replica
  if(length(technicalReplicas) > 1){
    biorepliaggregated <- aggregate(Intensity~Feature+Proteins+Condition+BioReplicate+Run, data = data2matrix, FUN = sum)
    biorepliaggregated <- aggregate(Intensity~Feature+Proteins+Condition+BioReplicate, data = biorepliaggregated, FUN = mean)
  }else{
    biorepliaggregated <- aggregate(Intensity~Feature+Proteins+Condition+BioReplicate, data = data2matrix, FUN = sum)
  }
  # Now let's sum up based on conditions
  biorepliaggregated <- aggregate(Intensity~Feature+Proteins+Condition, data = biorepliaggregated, FUN = median)
  biorepliaggregated$Intensity <- log2(biorepliaggregated$Intensity)
  evidencekeyscleanDCASTconditions <- dcast(data=biorepliaggregated, Proteins+Feature~Condition, value.var = "Intensity")
  precordfConditions <- evidencekeyscleanDCASTconditions[,3:dim(evidencekeyscleanDCASTconditions)[2]]
  Mcond <- cor(precordfConditions, use = "pairwise.complete.obs")
  
  theCondCorDis <- artms_plotCorrelationDistribution(Mcond)
  
  cat("\t~~> Conditions correlation matrix\n")
  matrixCorrelationCond <- gsub("evidence.txt", "correlationMatrixConditions.pdf", evidence_file)
  matrixCorrelationCond <- paste0("plot.",prot_exp,".",matrixCorrelationCond)
  matrixCorrelationCond <- paste0(output_dir,"/",matrixCorrelationCond)
  pdf(matrixCorrelationCond)
  corrplot(Mcond, 
           method="square", 
           addCoef.col = "white", 
           number.cex=0.6, 
           tl.col = "black", 
           title = "Matrix Correlation based on protein Intensities")
  pheatmap(Mcond, 
           cluster_rows = T,
           cluster_cols = T,
           cellheight = 10, 
           cellwidth=25, 
           main = "Clustering Conditions",
           fontsize=6, 
           fontsize_row=8, 
           fontsize_col=12, 
           border_color='black',
           fontfamily="Helvetica",
           treeheight_row = F, 
           treeheight_col = F,
           color = color.palette
  )
  print(theCondCorDis)
  garbage <- dev.off()
  
  # cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  # CONTRAST FILE REPRODUCIBILITY PLOTS
  # If a contrast file is provided it means that we want to generate reproducibility plots
  
  normalizeEqualizeMedians <- function(sedata){
    
    sedata <- biorepliaggregated
    # Get the list of Bioreplicates
    nbiorep <- unique(sedata$BioReplicate)
    
    # Choose the median that you want to use to normalize the data...
    # GLOBAL median
    medianChosen <- median(sedata$Intensity)
    
    ## Normalized by equilizing the Medians
    for (i in 1:length(nbiorep)) {
      medianBiorep <- median(sedata[sedata$BioReplicate == nbiorep[i], "Intensity"])
      # cat("\t\t\t",nbiorep[i], " -> ", medianBiorep,"\n")
      sedata[sedata$BioReplicate == nbiorep[i], "Intensity"] <- sedata[sedata$BioReplicate == nbiorep[i], "Intensity"]- medianBiorep + medianChosen
    }
    
    # Just check that truly happened...
    if(median(sedata$Intensity) != medianChosen){
      stop("\n\nOH NO. PROBLEMS WITH THE NORMALIZATION!!!\n\n\n")
    }else{
      cat("\tData Normalized on equilized medians\n\n")
      return(sedata)
    }
  }

  # DETAILS ABOUT THE SEARCH
  if (prot_exp == "apms" | prot_exp == "ab") {
    ekselect <- evidencekeysclean[c('Feature', 'Proteins', 'Intensity', 'Condition', 'BioReplicate', 'Run')]
    # Aggregate the technical replicas
    ekselecta <- aggregate(Intensity~Proteins+Condition+BioReplicate+Run, data=ekselect, FUN = sum)
    ekselectaBioreplica <- aggregate(Intensity~Proteins+Condition+BioReplicate+Run, data=ekselecta, FUN = sum) 
    # Select Unique number of proteins per Condition
    ac <- ekselecta[c('Proteins','Condition')]
    b <- unique(ac)
    # Select Unique number of proteins per Biological Replica
    bc <- ekselecta[c('Proteins','Condition','BioReplicate')]
    bb <- unique(bc)
    # Select unique number of proteins in Technical Replicates
    cc <- ekselecta[c('Proteins','Condition','BioReplicate', 'Run')]
    cc$TR <- paste0(ekselecta$BioReplicate,"_",ekselecta$Run)
    ccc <- cc[c('Proteins','Condition','TR')]
    cd <- unique(ccc)
  }else if (prot_exp == "ub") {
    ekselectall <- evidencekeysclean[c('Feature', 'Modified.sequence', 'Proteins', 'Intensity', 'Condition', 'BioReplicate', 'Run')]
    ekselectgly <- ekselectall[grep("(gl)", ekselectall$Modified.sequence),]
    # Select the REDUNDANT peptides with the meanimum intensity
    ekselectaBioreplica <- aggregate(Intensity~Proteins+Condition+BioReplicate+Run, data=ekselectgly, FUN = sum)
    # Select Unique number of proteins per Condition
    ac <- ekselectaBioreplica[c('Proteins','Condition')]
    b <- unique(ac)
    # Select Unique number of proteins per Biological Replica
    bc <- ekselectaBioreplica[c('Proteins','Condition','BioReplicate')]
    bb <- unique(bc)
    # Select unique number of proteins in Technical Replicates
    cc <- ekselectaBioreplica[c('Proteins','Condition','BioReplicate', 'Run')]
    cc$TR <- paste0(ekselectaBioreplica$BioReplicate,"_",ekselectaBioreplica$Run)
    ccc <- cc[c('Proteins','Condition','TR')]
    cd <- unique(ccc)
    
    # Check the total number of modified and non-modified residues to plot
    evidencekeys$MODIFICATION <- ifelse(grepl("(gl)",evidencekeys$Modified.sequence),"ub","other")
  }else if (prot_exp == "ph") {
    ekselectall <- evidencekeysclean[c('Feature', 'Modified.sequence', 'Proteins', 'Intensity', 'Condition', 'BioReplicate', 'Run')]
    ekselectgly <- ekselectall[grep("(ph)", ekselectall$Modified.sequence),]
    # SUM UP all the peptide intensities for all the proteins: only one intensity value per protein and biological data
    ekselectaBioreplica <- aggregate(Intensity~Proteins+Condition+BioReplicate+Run, data=ekselectgly, FUN = sum) 
    # Select Unique number of proteins per Condition
    ac <- ekselectaBioreplica[c('Proteins','Condition')]
    b <- unique(ac)
    # Select Unique number of proteins per Biological Replica
    bc <- ekselectaBioreplica[c('Proteins','Condition','BioReplicate')]
    bb <- unique(bc)
    # Select unique number of proteins in Technical Replicates
    cc <- ekselectaBioreplica[c('Proteins','Condition','BioReplicate', 'Run')]
    cc$TR <- paste0(ekselectaBioreplica$BioReplicate,"_",ekselectaBioreplica$Run)
    ccc <- cc[c('Proteins','Condition','TR')]
    cd <- unique(ccc)
    
    # Check the total number of modified and non-modified residues to plot
    evidencekeys$MODIFICATION <- ifelse(grepl("(ph)",evidencekeys$Modified.sequence),"ph","other")
    
  }else {
    stop("!-!-!-!- The proteomics experiment is not recognized (prot_exp should be apms, ab, ub, or ph)\n")
  }
  
  cat("\t~~>",prot_exp,"processed\n")
  
  cat("\t~~> Printing out the plots\n")
  reproName <- gsub("evidence.txt", "intensityStats.pdf", evidence_file)
  reproName <- paste0("plot.",prot_exp,".",reproName)
  reproName <- paste0(output_dir,"/",reproName)
  
  pdf(reproName)
  
  #QC: SUM of intensities per biological replica (peptides vs contaminant)
  ggplot(evisummary, aes(x=BioReplicate, y=Intensity, fill=Contaminant)) +
    geom_bar(stat="identity", position=position_dodge(width = 0.7), width=0.7) +
    theme_minimal() +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.title = element_blank()) +
    ggtitle("QC: Total Sum of Intensities in BioReplicates") +
    scale_fill_brewer(palette="Paired")
  
  #QC: SUM of intensities per condition (peptides vs contaminant)
  ggplot(evisummary, aes(x=Condition, y=Intensity, fill=Contaminant)) +
    geom_bar(stat="identity", position=position_dodge(width = 0.7), width=0.7) +
    theme_minimal() +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.title = element_blank()) +
    ggtitle("QC: Total Sum of Intensities in Conditions") +
    scale_fill_brewer(palette="Paired")
  
  #QC: TOTAL COUNT OF PEPTIDES IN EACH BIOLOGICAL REPLICA
  ggplot(evigeneral, aes(x = BioReplicate, fill = Contaminant)) + 
    geom_bar(stat = "count", position=position_dodge(width = 0.7), width=0.7) + 
    theme_minimal() +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.title = element_blank()) +
    ggtitle("QC: Total Peptide Counts in BioReplicates")
  
  #QC: TOTAL COUNT OF PEPTIDES IN EACH BIOLOGICAL REPLICA
  ggplot(evigeneral, aes(x = Condition, fill = Contaminant)) + 
    geom_bar(stat = "count", position=position_dodge(width = 0.7), width=0.7) + 
    theme_minimal() +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.title = element_blank()) +
    ggtitle("QC: Peptide Counts in Conditions")
  
  
  
  #QC: MAX INTENSITY IN REDUNDANT PEPTIDES, AND AGGREGATE FOR EACH PROTEIN THE SUM OF INTENSITY
  ekselectaBioreplica2plot <- ekselectaBioreplica
  ekselectaBioreplica2plot$Intensity <- log2(ekselectaBioreplica2plot$Intensity)
  ggplot(ekselectaBioreplica2plot, aes(x=as.factor(BioReplicate), y=Intensity, fill=Intensity))+
    geom_boxplot(aes(fill = BioReplicate))+
    # scale_y_log10() + coord_trans(y = "log10") +
    theme_linedraw() +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position = "none") +
    labs(x="BioReplicate", y="log2(Intensity)") + 
    ggtitle("Protein Intensity in BioReplicates\nExcluding contaminants. Max intensity of TR")
  
  ggplot(ekselectaBioreplica2plot, aes(x=Condition, y=Intensity, fill=Intensity))+
    geom_boxplot(aes(fill = Condition))+
    # scale_y_log10() + coord_trans(y = "log10") +
    theme_linedraw() +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position = "none") +
    labs(x="Condition", y="log2(Intensity)") + 
    ggtitle("Protein Intensity in Conditions\nExcluding contaminants. Max intensity of TR")
  
  ggplot(ekselectaBioreplica) +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
    geom_bar(aes(BioReplicate, Intensity), position = "dodge", stat = "summary", fun.y = "mean", fill="black", colour = "orange") +
    ggtitle("Total Intensity in Biological Replicas\nExcluding contaminants. Max intensity of TR")  
  
  ggplot(ekselectaBioreplica) +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
    geom_bar(aes(Condition, Intensity), position = "dodge", stat = "summary", fun.y = "mean", fill="black", colour = "green") +
    ggtitle("Total Intensity in Conditions\nExcluding contaminants. Max intensity of TR")  
  
  ggplot(cd, aes(x = TR, fill = Condition)) + 
    geom_bar(stat = "count") + 
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.9), legend.position = "none") +
    geom_text(stat='count', aes(label=..count..), vjust=-0.5, size = 2.7) +
    ggtitle("Unique Features in Technical Replicas")    
  
  ggplot(bb, aes(x = BioReplicate, fill = Condition)) + 
    geom_bar(stat = "count") + 
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.position = "none") +
    geom_text(stat='count', aes(label=..count..), vjust=-0.5, size = 2.7) +
    ggtitle("Unique Features in Biological Replicas")
  
  ggplot(b, aes(x = Condition, fill = Condition)) + 
    geom_bar(stat = "count") + 
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.position = "none") +
    geom_text(stat='count', aes(label=..count..), vjust=-0.5, size = 2.7) +
    ggtitle("Unique Features in Condition")    
  garbage <- dev.off()
  
  if( prot_exp == "ph" | prot_exp == "ub" | prot_exp == "ptmph"){
    cat("\t~~> Printing",prot_exp,"count\n")
    modName <- gsub("evidence.txt", "ptmStats.pdf", evidence_file)
    modName <- paste0("plot.",prot_exp,".",modName)
    modName <- paste0(output_dir,"/",modName)
    
    x <- ggplot(evidencekeys, aes(x = BioReplicate, fill = MODIFICATION))
    x <- x + geom_bar(stat = "count", position=position_dodge(width = 0.7), width=0.7)
    x <- x + theme_minimal()
    x <- x + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.title = element_blank())
    x <- x + ggtitle("Peptide Count in Biological Replicas")
    
    y <- ggplot(evidencekeys, aes(x = Condition, fill = MODIFICATION))
    y <- y + geom_bar(stat = "count", position=position_dodge(width = 0.7), width=0.7)
    y <- y + theme_minimal()
    y <- y + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.title = element_blank())
    y <- y + ggtitle("Peptide Count in Conditions")
    
    u <- ggplot(evidencekeys, aes(x = BioReplicate, y=Intensity, fill = MODIFICATION))
    u <- u + geom_bar(stat = "identity", position=position_dodge(width = 0.7), width=0.7)
    u <- u + theme_minimal()
    u <- u + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.title = element_blank())
    u <- u + ggtitle("Total Peptide Intensity in Biological Replicas")
    u <- u + scale_fill_brewer(palette="Paired")
    
    z <- ggplot(evidencekeys, aes(x = Condition, y=Intensity, fill = MODIFICATION))
    z <- z + geom_bar(stat = "identity", position=position_dodge(width = 0.7), width=0.7)
    z <- z + theme_minimal()
    z <- z + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.title = element_blank())
    z <- z + ggtitle("Total Peptide Intensity in Conditions")
    z <- z + scale_fill_brewer(palette="Paired")
    
    #QC: Total Peptides
    pdf(modName)
      print(x)
      print(y)
      print(u)
      print(z)
    garbage <- dev.off()
  }
  
  cat("\t~~> Done!\n\n")
}



## SITE

# setwd('~/Box Sync/projects/krogan/judd/Roche/HBV/JFH102_virus/')
# evidence_file <- 'roche-jfh102-papd57-yeshbv-evidence.txt'
# keys_file <- 'roche-jfh102-papd57-yeshbv-keys.txt'
# prot_exp <- 'apms'
# output_dir <- 'testing123'
# fractions <- 'nofra'
# contrast <- 'roche-jfh102-papd57-yeshbv-contrast.txt'
# specie <- 'human'

# setwd('~/Box Sync/projects/sanfordBurnham/Deshpande/data/')
# evidence_file <- 'deshpande-evidence.txt'
# keys_file <- 'deshpande-keys.txt'
# prot_exp <- 'ab'

# setwd('~/Box Sync/20171020_CH1/')
# evidence_file <- 'evidence.txt'
# keys_file <- 'keys.txt'
# prot_exp <- 'ph'


# ZIKA
# setwd('~/Box Sync/projects/krogan/zika/mist_analysis/data/')
# evidence_file <- '20180313-zika-krystal-evidence.txt'
# keys_file <- '20180313-zika-krystal-keys.txt'
# prot_exp <- 'apms'
# output_dir <- 'whatever'
# fractions <- 'nofra'
# contrast <- 'nocon'


# # Siah2 Zeev
# setwd('~/Box Sync/danielleSwaney/ZeevRonaiLab/DS1/ab/')
# evidence_file <- 'ronai_siah2_ab-evidence.txt'
# keys_file <- 'ronai_siah2_ab-keys.txt'
# prot_exp <- 'ab'
# output_dir <- 'whatever'
# fractions <- 'nofra'
# contrast <- 'nocon'
# specie <- 'mouse'
# 
# setwd('~/Box Sync/danielleSwaney/ZeevRonaiLab/DS1/ub/')
# evidence_file <- 'ronai_siah2_ub-evidence.txt'
# keys_file <- 'ronai_siah2_ub-keys.txt'
# prot_exp <- 'ub'
# output_dir <- 'whatever'
# fractions <- 'nofra'
# contrast <- 'nocon'