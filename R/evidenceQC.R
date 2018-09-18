# ------------------------------------------------------------------------------
#' @title Quality Control analysis of the MaxQuant evidence file
#' 
#' @description Quality Control analysis of the MaxQuant evidence file
#' @param evidence_file (char) The evidence file path and name
#' @param keys_file (char) The keys file path and name
#' @param output_name (char) prefix output name (no extension). 
#' Default: "qcPlots_evidence"
#' @param prot_exp (char) Proteomics experiment. 4 options available:
#' - `APMS`: affinity purification mass spectrometry
#' - `AB`: protein abundance
#' - `PH`: protein phosphorylation
#' - `UB`: protein ubiquitination (aka ubiquitylation)
#' @param fractions (binary) Is a fractionated experiment?
#' - 1 yes
#' - 0 no (default)
#' @return Quality control files and plots
#' @keywords QC, quality, control, evidence
#' @examples{
#' artms_evidenceQC(evidence_file = artms_data_ph_evidence, 
#'                  keys_file = artms_data_ph_keys, 
#'                  output_name = "qcPlots_evidence",
#'                  prot_exp = "PH")
#' }
#' @export
artms_evidenceQC <- function(evidence_file,
                             keys_file, 
                             prot_exp,
                             fractions = 0, 
                             output_name = "qcPlots_evidence"){

  # evidence_file <- Sys.glob(evidence_file)
  # keys_file <- Sys.glob(keys_file)
  
  prot_exp <- toupper(prot_exp)

  if(any(!prot_exp %in% c('AB','PH','UB','APMS'))){
    cat("\nERROR!!!\nTHE prot_exp ARGUMENT IS NOT CORRECT.\n")
    cat("IT MUST BE ONE OF THE FOLLOWINGS:\n\t- AB\n\t- PH\n\t- UB\n\t- APMS\n")
    stop("PLEASE, PROVIDE A CORRECT prot_exp ARGUMENT")
  }
  
  if(fractions){
    # Check that the keys file is correct
    keys <- .artms_checkIfFile(keys_file)
    keys <- .artms_checkRawFileColumnName(keys)
    if(any(!'FractionKey' %in% colnames(keys))){
      cat('\nERROR!!! fractions WAS ACTIVATED BUT FractionKey COLUMN NOT FOUND IN THE KEYS FILE\n')
      stop('Please, try again once revised\n\n')
    }
  }
  
  cat("\nQUALITY CONTROL ------------\n")
  cat(">> LOADING THE EVIDENCE FILE\n")
  cat("(it should take some time due to the usual large size of evidence files)\n")
  
  evidence <- .artms_checkIfFile(evidence_file)
  evidence <- .artms_checkRawFileColumnName(evidence)
  
  keys <- .artms_checkIfFile(keys_file)
  keys <- .artms_checkRawFileColumnName(keys)

  # EVIDENCE:
  evidencekeys <- artms_mergeMaxQDataWithKeys(evidence, keys)
  
  ekselecta <- aggregate(Intensity~Proteins+Condition+BioReplicate+Run, data=evidencekeys, FUN = sum)
  ekselectaBioreplica <- aggregate(Intensity~Proteins+Condition+BioReplicate, data=ekselecta, FUN = sum) 
  
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
  
  cat(">> GENERATING THE INTENSITY DISTRIBUTION PLOTS\n")
  intDistribution <- paste0(output_name,".qcplot.IntensityDistributions.pdf")
  
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
    evidencekeys <- artms_changeColumnName(evidencekeys,"Leading.Proteins","Leading.proteins")
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
  
  # AGGREGATE ALL THE INTENSITIES PER PROTEIN, summing everything up
  evigeneral$Intensity <- as.numeric(evigeneral$Intensity)
  
  evisummary <- aggregate(Intensity~Condition+BioReplicate+Contaminant+Run, data=evigeneral, FUN = sum)
  
  # CLEANING THE EVIDENCE OF CONTAMINANTS
  evidencekeysclean <- artms_filterMaxqData(evidencekeys)
  
  cat(">> GENERATING THE REPRODUCIBILITY PLOTS (warning: it will take some time)\n")
  seqReproName <- paste0(output_name,".qcplot.basicReproducibility.pdf")
  
  if (prot_exp == "UB") {
    evidencekeysclean <- evidencekeysclean[c('Feature', 'Modified.sequence', 'Proteins', 'Intensity', 'Condition', 'BioReplicate', 'Run')]
    evidencekeysclean <- evidencekeysclean[grep("(gl)", evidencekeysclean$Modified.sequence),]  
  }else if (prot_exp == "PH") {
    evidencekeysclean <- evidencekeysclean[c('Feature', 'Modified.sequence', 'Proteins', 'Intensity', 'Condition', 'BioReplicate', 'Run')]
    evidencekeysclean <- evidencekeysclean[grep("(ph)", evidencekeysclean$Modified.sequence),]
  }
  
  pdf(seqReproName)
    .artms_plotReproducibilityEvidence(evidencekeysclean)
  garbage <- dev.off()
  
  # Create matrix of reproducibility TECHNICAL REPLICAS
  data2matrix <- evidencekeysclean
  # Make sure that the Intensity is numeric
  data2matrix$Intensity <- as.numeric(data2matrix$Intensity)
  
  # Check the number of TECHNICAL REPLICAS by checking the first technical replica
  technicalReplicas <- unique(data2matrix$Run[which(data2matrix$BioReplicate == data2matrix$BioReplicate[1])])
  palette.breaks <- seq(1, 3, 0.1)
  color.palette  <- colorRampPalette(c("white","steelblue"))(length(palette.breaks))
  
  cat(">> GENERATING CORRELATION MATRICES\n")
  if(length(technicalReplicas) > 1){
    # First aggregate at the protein level by summing up everything
    biorepliaggregated <- aggregate(Intensity~Feature+Proteins+Condition+BioReplicate+Run, data = data2matrix, FUN = sum)
    biorepliaggregated$Intensity <- log2(biorepliaggregated$Intensity)
    evidencekeyscleanDCASTbioreplicas <- data.table::dcast(data=biorepliaggregated, Proteins+Feature~BioReplicate+Run, value.var = "Intensity")
    precordfBioreplicas <- evidencekeyscleanDCASTbioreplicas[,3:dim(evidencekeyscleanDCASTbioreplicas)[2]]
    Mtechnicalrep <- cor(precordfBioreplicas, use = "pairwise.complete.obs")
    
    theTechCorDis <- .artms_plotCorrelationDistribution(Mtechnicalrep)
    
    # And now for clustering
    cat("--- By Technical replicates\n")
    matrixCorrelationBioreplicas <- paste0(output_name, ".qcplot.correlationMatrixTR.pdf")
    
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
    cat("--- NO Technical Replicates detected\n")
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
  evidencekeyscleanDCASTbioreplicas <- data.table::dcast(data=biorepliaggregated, Proteins+Feature~BioReplicate, value.var = "Intensity")
  precordfBioreplicas <- evidencekeyscleanDCASTbioreplicas[,3:dim(evidencekeyscleanDCASTbioreplicas)[2]]
  Mbioreplicas <- cor(precordfBioreplicas, use = "pairwise.complete.obs")
  
  theBiorCorDis <- .artms_plotCorrelationDistribution(Mbioreplicas)
  
  cat("--- By Biological replicates\n")
  matrixCorrelationBioreplicas <- paste0(output_name,".qcplot.correlationMatrixBR.pdf")
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
  evidencekeyscleanDCASTconditions <- data.table::dcast(data=biorepliaggregated, Proteins+Feature~Condition, value.var = "Intensity")
  precordfConditions <- evidencekeyscleanDCASTconditions[,3:dim(evidencekeyscleanDCASTconditions)[2]]
  Mcond <- cor(precordfConditions, use = "pairwise.complete.obs")
  
  theCondCorDis <- .artms_plotCorrelationDistribution(Mcond)
  
  cat("--- By Conditions\n")
  matrixCorrelationCond <- paste0(output_name, ".qcplot.correlationMatrixConditions.pdf")
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
  
  # DETAILS
  cat(">> GENERATING INTENSITY STATS PLOTS\n")
  if (prot_exp == "APMS" | prot_exp == "AB") {
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
  }else if (prot_exp == "UB") {
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
  }else if (prot_exp == "PH") {
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
    stop("PROTEOMICS EXPERIMENT NOT RECOGNIZED (prot_exp SHOULD BE APMS, AB, UB, or PH)\n")
  }
  
  cat("--- ",prot_exp," PROCESSED\n")
  reproName <- paste0(output_name, ".qcplot.intensityStats.pdf")

  #QC: SUM of intensities per biological replica (peptides vs contaminant)
  pisa <- ggplot(evisummary, aes(x=BioReplicate, y=Intensity, fill=Contaminant)) +
    geom_bar(stat="identity", position=position_dodge(width = 0.7), width=0.7) +
    theme_minimal() +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.title = element_blank()) +
    ggtitle("QC: Total Sum of Intensities in BioReplicates") +
    scale_fill_brewer(palette="Paired")
  
  #QC: SUM of intensities per condition (peptides vs contaminant)
  pisb <- ggplot(evisummary, aes(x=Condition, y=Intensity, fill=Contaminant)) +
    geom_bar(stat="identity", position=position_dodge(width = 0.7), width=0.7) +
    theme_minimal() +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.title = element_blank()) +
    ggtitle("QC: Total Sum of Intensities in Conditions") +
    scale_fill_brewer(palette="Paired")
  
  #QC: TOTAL COUNT OF PEPTIDES IN EACH BIOLOGICAL REPLICA
  pisc <- ggplot(evigeneral, aes(x = BioReplicate, fill = Contaminant)) + 
    geom_bar(stat = "count", position=position_dodge(width = 0.7), width=0.7) + 
    theme_minimal() +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.title = element_blank()) +
    ggtitle("QC: Total Peptide Counts in BioReplicates")
  
  #QC: TOTAL COUNT OF PEPTIDES IN EACH BIOLOGICAL REPLICA
  pisd <- ggplot(evigeneral, aes(x = Condition, fill = Contaminant)) + 
    geom_bar(stat = "count", position=position_dodge(width = 0.7), width=0.7) + 
    theme_minimal() +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.title = element_blank()) +
    ggtitle("QC: Peptide Counts in Conditions")
  
  #QC: MAX INTENSITY IN REDUNDANT PEPTIDES, AND AGGREGATE FOR EACH PROTEIN THE SUM OF INTENSITY
  ekselectaBioreplica2plot <- ekselectaBioreplica
  ekselectaBioreplica2plot$Intensity <- log2(ekselectaBioreplica2plot$Intensity)
  pise <- ggplot(ekselectaBioreplica2plot, aes(x=as.factor(BioReplicate), y=Intensity, fill=Intensity))+
    geom_boxplot(aes(fill = BioReplicate))+
    # scale_y_log10() + coord_trans(y = "log10") +
    theme_linedraw() +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position = "none") +
    labs(x="BioReplicate", y="log2(Intensity)") + 
    ggtitle("Protein Intensity in BioReplicates\nExcluding contaminants. Max intensity of TR")
  
  pisf <- ggplot(ekselectaBioreplica2plot, aes(x=Condition, y=Intensity, fill=Intensity))+
    geom_boxplot(aes(fill = Condition))+
    # scale_y_log10() + coord_trans(y = "log10") +
    theme_linedraw() +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position = "none") +
    labs(x="Condition", y="log2(Intensity)") + 
    ggtitle("Protein Intensity in Conditions\nExcluding contaminants. Max intensity of TR")
  
  pisg <- ggplot(ekselectaBioreplica) +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
    geom_bar(aes(BioReplicate, Intensity), position = "dodge", stat = "summary", fun.y = "mean", fill="black", colour = "orange") +
    ggtitle("Total Intensity in Biological Replicas\nExcluding contaminants. Max intensity of TR")  
  
  pish <- ggplot(ekselectaBioreplica) +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
    geom_bar(aes(Condition, Intensity), position = "dodge", stat = "summary", fun.y = "mean", fill="black", colour = "green") +
    ggtitle("Total Intensity in Conditions\nExcluding contaminants. Max intensity of TR")  
  
  pisi <- ggplot(cd, aes(x = TR, fill = Condition)) + 
    geom_bar(stat = "count") + 
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.9), legend.position = "none") +
    geom_text(stat='count', aes(label=..count..), vjust=-0.5, size = 2.7) +
    ggtitle("Unique Features in Technical Replicas")    
  
  pisj <- ggplot(bb, aes(x = BioReplicate, fill = Condition)) + 
    geom_bar(stat = "count") + 
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.position = "none") +
    geom_text(stat='count', aes(label=..count..), vjust=-0.5, size = 2.7) +
    ggtitle("Unique Features in Biological Replicas")
  
  pisk <- ggplot(b, aes(x = Condition, fill = Condition)) + 
    geom_bar(stat = "count") + 
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.position = "none") +
    geom_text(stat='count', aes(label=..count..), vjust=-0.5, size = 2.7) +
    ggtitle("Unique Features in Condition")
  
  pdf(reproName)
    print(pisa)
    print(pisb)
    print(pisc)
    print(pisd)
    print(pise)
    print(pisf)
    print(pisg)
    print(pish)
    print(pisi)
    print(pisj)
    print(pisk)
  garbage <- dev.off()
  
  if( prot_exp == "PH" | prot_exp == "UB"){
    cat(">> GENERATING PTM ",prot_exp," STATS\n")
    modName <- paste0(output_name, "qcplot.ptmStats.pdf")
    
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
  
  cat("\n>>QUALITY CONTROL ANALYSIS COMPLETED!\n")
}

