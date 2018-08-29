# ------------------------------------------------------------------------------
#' @title Analysis of the Relative Quantification
#' 
#' @description Analysis of the Relative Quantifications obtained by MSstats.
#' It includes:
#' 
#' - Annotations
#' @param log2fc_file MSstats results file 
#' @param modelqc_file MSstats modelqc file
#' @param specie Specie (human, mouse)
#' @param enrich Performed enrichment analysis?
#' @param output_dir results folder name
#' @param isFluomics Is from the fluomics project?
#' @param isPtm is a ptm quantification?
#' @param isBackground background gene set
#' @param mnbr minimal number of biological replicates for imputation
#' @param threshold log2fc cutoff for enrichment analysis
#' @param ipval pvalue cutoff for enrichment analysis
#' @param pathogen is there a pathogen in the dataset as well?
#' @return summary of quantifications, including annotations, enrichments, etc
#' @keywords analysis, quantifications
#' artms_analysisQuantifications()
#' @export
artms_analysisQuantifications <- function(log2fc_file, 
                                          modelqc_file, 
                                          specie, 
                                          rm_contaminant, 
                                          enrich, 
                                          output_dir, 
                                          isFluomics, 
                                          isPtm, 
                                          isBackground, 
                                          mnbr, 
                                          threshold, 
                                          ipval, 
                                          pathogen){
  
  log2fc_file = "ab-testing-new-results.txt"
  modelqc_file = "ab-testing-new-results_ModelQC.txt"
  specie = "human"
  enrich = "yes"
  output_dir = "resultsTesting"
  isFluomics = "yes"
  isPtm = "noptmsites"
  isBackground = "nobackground"
  mnbr = 2
  threshold = 1
  ipval = "pvalue"
  pathogen = "nopathogen"
  
  # source('~/github/kroganlab/enrichment/enrichProfiler.R')
  # source('~/github/kroganlab/djmSource/myLibrary/generalFunctions.R')
  
  # It required all these installations:
  # source("https://bioconductor.org/biocLite.R")
  # biocLite("MSstats")
  # install.packages('openxlsx')
  # install.packages('pheatmap')
  # install.packages('tidyr')
  # biocLite("msa")
  # install.packages('FactoMineR')
  # install.packages('devtools')
  # install.packages('factoextra')
  # install.packages('corrplot')
  # install.packages('PerformanceAnalytics')
  
  # suppressMessages(library(reshape2))
  # suppressMessages(library(openxlsx))
  # suppressMessages(library(pheatmap))
  # suppressMessages(library(ggplot2))
  # suppressMessages(library(ggrepel))
  # suppressMessages(library(tidyr))
  # suppressMessages(library(dplyr))
  
  cat(">> ANALYSIS OF QUANTIFICATIONS\n")
  
  if(pathogen == "nopathogen"){
    cat("No Pathogen extra in these samples (choose this for Influenza)\n")
  }else if(pathogen == "tb"){ # This should not work
    cat("PATHOGEN IN SAMPLES: TB\n")
    pathogen.ids <- read.delim('~/Box Sync/db/uniprot/uniprot-tr-myctb_tuberculosis_ATCC35801_TMC10-onlyEntryID.fasta', header = F, sep = "\t", quote = "", stringsAsFactors = F) # pathogen.ids$Entry, "TB",
    names(pathogen.ids) <- c('Entry')
  }else if(pathogen == "lpn"){
    cat("PATHOGEN IN SAMPLES: LEGIONELLA PNEUMOPHILA\n")
    pathogen.ids <- read.delim('~/Box Sync/db/uniprot/uniprot-legionella-proteome_UP000000609.txt', header = T, sep = "\t", quote = "", stringsAsFactors = F) # pathogen.ids$Entry, "Lpn",
  }else{
    stop("\n\nThis pathogen is not supported yet\n\n")
  }
  
  output_dir <- paste0(output_dir,"_",ipval)
  
  # create output directory if it doesn't exist
  if(!dir.exists(output_dir)){
    dir.create(output_dir, recursive = T)
  }
  
  # LOADING ABUNDANCE
  cat("---READING IN modelqc FILE\n")
  dfmq <- read.delim(modelqc_file, header = T, sep = "\t", stringsAsFactors = F)
  #Removing the empty protein names
  if(any(dfmq$PROTEIN == "")){ dfmq <- dfmq[-which(dfmq$PROTEIN == ""),]}
  dfmq$PROTEIN <- gsub("(sp\\|)(.*)(\\|.*)", "\\2", dfmq$PROTEIN )
  dfmq$PROTEIN <- gsub("(.*)(\\|.*)", "\\1", dfmq$PROTEIN )
  # Select only POSITIVE values (they all should be positive values, 
  # but in one case we found negative ones)
  dfmq <- dfmq[which(dfmq$ABUNDANCE > 12 & dfmq$ABUNDANCE < 35),]
  
  # First, let's take the conditions, which will be used later in several places
  conditions <- unique(dfmq$GROUP_ORIGINAL)
  numberConditions <- length(conditions)
  

  # KEY STEP: GETTING THE BACKGROUND GENE LIST
  if(isBackground == "nobackground"){
    # If not list of background genes is provided, 
    # then extract them from the modelqc file
    if (isPtm == "noptmsites"){
      dfmq2Genes <- artms_annotationUniprot(dfmq, 'PROTEIN', specie)
      numberTotalGenes <- length(unique(dfmq2Genes$Gene))
      cat(">> TOTAL NUMBER OF GENES/PROTEINS: ",numberTotalGenes,"\n")
      listOfGenes <- unique(dfmq2Genes$Gene)
    }else if( grepl("yesptm",isPtm) ){
      dfmq2Genes <- dfmq[c('PROTEIN','GROUP_ORIGINAL')] # If you want to apply some sort of filter, do it here
      names(dfmq2Genes)[grep('PROTEIN',names(dfmq2Genes))] <- 'Protein'
      # Removing party sites
      dfmq2Genes <- dfmq2Genes[grep(",",dfmq2Genes$Protein, invert=T),]
      # And now be very careful with the Fluomics labeling, since they have an extra _ that it is not follow by the site
      cat("---Warning! if you have protein_sites id with more than one '_' is going to be a problem\n\n")
      dfmq2Genes$Protein <- ifelse(
        grepl("_H1N1|_H3N2|_H5N1", dfmq2Genes$Protein), gsub("^(\\S+?_H[1,3,5]N[1,2])_.*", "\\1", dfmq2Genes$Protein, perl = T) , 
        gsub("^(\\S+?)_.*", "\\1", dfmq2Genes$Protein, perl = T)
      )
      dfmq2Genes <- unique(dfmq2Genes)
      dfmq2Genes <- artms_annotationUniprot(dfmq2Genes, 'Protein', specie)
      numberTotalGenes <- length(unique(dfmq2Genes$Gene))
      cat(">> TOTAL NUMBER OF GENES/PROTEINS: ",numberTotalGenes,"\n\n")
      if(numberTotalGenes == 0){
        stop("\nSOMETHING WRONG WITH THE IDs. Check the source code\n")
      }
      listOfGenes <- unique(dfmq2Genes$Gene)
    }
  }else{
    # No matter what list is provided, it must come with a "Gene" column
    backgroundList <- read.delim(isBackground, header = T, sep = "\t", quote = "", stringsAsFactors = F)
    listOfGenes <- unique(backgroundList$Gene)
  }
  
  backgroundNumber <- length(listOfGenes)
  
  ##############################################################################
  # LOG2FC
  dflog2fcraw <- read.delim(log2fc_file, header = T, sep = "\t", stringsAsFactors = F)
  if(any(dflog2fcraw$Protein == "")){ dflog2fcraw <- dflog2fcraw[-which(dflog2fcraw$Protein == ""),]}
  dflog2fcraw$Protein <- gsub("(sp\\|)(.*)(\\|.*)", "\\2", dflog2fcraw$Protein )
  dflog2fcraw$Protein <- gsub("(.*)(\\|.*)", "\\1", dflog2fcraw$Protein )
  
  # # Filtering conditions:
  # if(isFluomics == "yesflu"){
  #   cat("WARNING! selecting MOCK for humans and mice in LOG2FC raw data\n")
  # 
  #   # WHEN A REFERENCE MOCK IS WISHED
  #   # if(specie == "human"){
  #   #   dflog2fcraw <- dflog2fcraw[(grepl("H[[:digit:]]N[[:digit:]]", dflog2fcraw$Label) & grepl("MOCK_03H", dflog2fcraw$Label) ),]
  #   # }else if (specie == "mouse"){
  #   #   dflog2fcraw <- dflog2fcraw[(grepl("H[[:digit:]]N[[:digit:]]", dflog2fcraw$Label) | grepl("MOCK_D04", dflog2fcraw$Label) ),]
  #   # }
  # 
  #   #Choosing the match mock
  #   flu_contrast <- c("H1N1_03H-MOCK_03H", "H1N1_06H-MOCK_06H", "H1N1_12H-MOCK_12H", "H1N1_18H-MOCK_18H", "H3N2_03H-MOCK_03H", "H3N2_06H-MOCK_06H", "H3N2_12H-MOCK_12H", "H3N2_18H-MOCK_18H", "H5N1_03H-MOCK_03H", "H5N1_06H-MOCK_06H", "H5N1_12H-MOCK_12H", "H5N1_18H-MOCK_18H")
  #   dflog2fcraw <- dflog2fcraw[which(dflog2fcraw$Label %in% flu_contrast),]
  # 
  #   cat("LOG2FC Data Filtered by specific FLU comparisons\n")
  # }
  
  # Let's get rid of outliers: log2fc larger than X (but we need to keep the "inf" values for imputation)
  dflog2fcfinites <- dflog2fcraw[is.finite(dflog2fcraw$log2FC),]
  cutofflog2fc <- 12
  filtermorethan10 <- length(dflog2fcfinites$log2FC[abs(dflog2fcfinites$log2FC) > cutofflog2fc])
  if( filtermorethan10 > 0){
    cat("\n\t Removing log2fc values larger than +/-",cutofflog2fc," ... ", filtermorethan10, "\n\n")
    dflog2fcfinites <- dflog2fcfinites[-which(abs(dflog2fcfinites$log2FC) > cutofflog2fc),]
  }else{
    cat("No log2fc values larger than 12, so moving on!\n")
  }
  
  # IMPUTING MISSING VALUES
  # When a value is completely missed in one of the conditions, 
  # the log2fc = Inf / -Inf. Here, we impute those values.
  # The imputation method works as follow. The assumption is that those
  # proteins are likely present as well in those conditions where are missed, but due to the
  # small sampling (usually 2 or 3 biological replicas) and other proteomics
  # related issue, those proteins didn't make it through the level of detection.
  # Therefore, a small intensity (sampled from the bottom 5%) will be assigned
  # to the protein/site in the missing condition, and the new log2fc is re-calculated
  # out of the MSstats box. Two issues are addressed in this way
  # 1. If a protein has been consistently identified in one of the conditions, it will stay
  # 2. But if the intensity value in those conditions was too low, then the log2fc will be also low
  
  # Select infinite values (i.e., log2fc missed for that)
  dflog2fcinfinites <- dflog2fcraw[is.infinite(dflog2fcraw$log2FC),]
  numberInfinites <- dim(dflog2fcinfinites)[1]
  
  # Control
  if( numberInfinites < 1){
    cat("Number of infinite values: \n")
    cat("\t\t Less than --> 1 infinite value\n")
    stop("\n\n\tQuality control stop: there are less than 1 missing values in this sample!\nCheck what's going on\n\n")
  }else {
    cat("Selecting infinitive values for imputation\n")
    cat("\t---> Number of infinite values",dim(dflog2fcinfinites)[1],"\n")
    # CONTROL CODE
    # The log2fc values are randomly calculated between two values at the bottom.
    # But just in case, next a plot will be generated to show that the log2fc values
    # for the imputed proteins/sites just slightly change if the whole thing is 
    # recalculated
    
    # # Loop to impute the values 100 times:
    # for(i in 1:100){
    #   here <- imputeMissingValues(dflog2fcinfinites, dfmq, i)
    #   names(here)[grep("log2FC", names(here))] <- paste0("log2fc",i)
    #   if(i == 1){
    #     tmp <- here
    #   }else{
    #     tmp <- merge(tmp,here, by = c('Protein','Label')) 
    #   }
    # }
    # 
    # # Select one of the comparisons:
    # onlyab <- tmp[which(tmp$Label == "A-B"),]
    # tmpsorted <- onlyab[order(onlyab$log2fc1),]
    # # Select top and bottom 10 and bind them
    # top10 <- head(tmpsorted, n = 10)
    # bottom10 <- tail(tmpsorted, n = 10)
    # tb10 <- rbind(top10, bottom10)
    # 
    # # Melt 2 PLOT:
    # melttop10 <- melt(tmpsorted, id.vars = c('Protein','Label'), variable.name = "serie", value.name = 'log2fc')
    # library(directlabels)
    # eicoEnzymes <- ggplot(melttop10, aes(serie, log2fc, group = Protein, colour = Protein)) +
    #   geom_point() + geom_line(alpha = 0.8) +
    #   scale_colour_discrete(guide = 'none') +
    #   scale_x_discrete(expand=c(0, 1)) +
    #   geom_dl(aes(label = Protein), method = list(dl.combine("first.points", "last.points"), cex = 0.9)) #library(directlabels)
    # eicoEnzymes <- eicoEnzymes + ggtitle("Imputing log2fc 100 times")
    # eicoEnzymes <- eicoEnzymes + ylim(16,-16)
    # eicoEnzymes <- eicoEnzymes + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position = "none")
    # pdf('testingtesting.pdf')
    # print(eicoEnzymes)
    # dev.off()
    
    imputedL2FCmelted <- imputeMissingValues(dflog2fcinfinites, dfmq)
    
    # Merge with the original log2fc values to impute...
    theImputedL2FC <- merge(dflog2fcinfinites, imputedL2FCmelted, by = c("Protein","Label"), all.x = T)
    theImputedL2FC$imputed <- "yes"
  }
  
  # Getting the data ready for merging
  dflog2fcfinites$imputed <- "no"
  dflog2fcfinites$iLog2FC <- dflog2fcfinites$log2FC
  
  # Choose the pvalue or adjusted pvalue as the iPvalue
  if(ipval == "pvalue"){
    dflog2fcfinites$iPvalue <- dflog2fcfinites$pvalue
  }else if(ipval == "adjpvalue"){
    dflog2fcfinites$iPvalue <- dflog2fcfinites$adj.pvalue
  }else{
    stop("\n\n\t------> wait a minute: did you choose pvalue or adjpvalue")
  }
  
  # Merging: NA values are thrown away at this point
  if( numberInfinites < 1){
    dflog2fc <- dflog2fcfinites
  } else {
    dflog2fc <- rbind(dflog2fcfinites, theImputedL2FC)
  }
  
  cat("Plotting distributions\n")
  
  plotDFdistColor <- ggplot(dflog2fc, aes(x = log2FC, fill = Label)) + 
    geom_histogram(bins = 100, alpha = .4, col="black") +
    labs(title="Distribution log2FC", x = "log2FC")
  
  plotDFdistAll <- ggplot(dflog2fc, aes(x = log2FC)) + 
    geom_histogram(bins = 100, alpha = .4, col="black") +
    labs(title="Distribution log2FC", x = "log2FC")
  
  plotDFdistiLog <- ggplot(dflog2fc, aes(x = iLog2FC)) + 
    geom_histogram(bins = 100, alpha = .4, col="black") +
    labs(title="Distribution ilog2FC (imputed + nonimputed", x = "iLog2FC")
  
  plotPvalues <- ggplot(dflog2fc[is.finite(dflog2fc$pvalue),], aes(x = pvalue)) +
    geom_histogram(bins = 50, alpha = .4, col = "black") +
    labs(title="Distribution p-values", x = "p-values")
  
  plotAdjustedPvalues <- ggplot(dflog2fc[-which(dflog2fc$adj.pvalue == 0),], aes(x = adj.pvalue)) + 
    geom_histogram(bins = 150, alpha = .4, col="black") +
    labs(title="Distribution adj.pvalues", x = "adj.values")
  
  plotAdjustedIpvalues <- ggplot(dflog2fc, aes(x = iPvalue)) + 
    geom_histogram(bins = 150, alpha = .4, col="black") +
    labs(title="Distribution imputed p-values", x = "iPvalues")
  
  
  # DISTRIBUTION PRINT OUTS
  distributionsOut <- gsub(".txt",".distributions.pdf",log2fc_file)
  distributionsOut <- paste0(output_dir,"/",distributionsOut)
  pdf(distributionsOut)
  plotDFdistColor
  plotDFdistAll
  plotDFdistiLog
  plotPvalues
  plotAdjustedPvalues
  plotAdjustedIpvalues
  
  if( numberInfinites > 0){
    hist(imputedL2FCmelted$iLog2FC, breaks = 100, main = paste0("Imputed Log2FC (all)\n  n = ",dim(imputedL2FCmelted)[1]), xlab = "log2fc")
    hist(theImputedL2FC$iLog2FC, breaks = 100, main = paste0("Imputed Log2FC merged\n n = ", dim(theImputedL2FC)[1] ), xlab = "log2fc")
  }
  hist(dflog2fcfinites$pvalue, breaks = 100, main = paste0("p-value distribution\n n = ",dim(dflog2fcfinites)[1]), xlab = "adj.pvalues")
  hist(dflog2fcfinites$adj.pvalue, breaks = 100, main = paste0("Adjusted p-values distribution\n n = ", dim(dflog2fcfinites)[1]), xlab = "adj.pvalues")
  hist(dflog2fcfinites$iLog2FC, breaks = 1000, main = paste0("Non-imputed Log2FC distribution\n n = ", dim(dflog2fcfinites)[1]), xlab = "log2FC")
  hist(dflog2fc$iPvalue, breaks = 100, main = paste0("(Imputed+NonImputed) adjusted pvalue distribution\n n = ",dim(dflog2fc)[1]), xlab = "adj.pvalues")
  hist(dflog2fc$iLog2FC, breaks = 1000, main = paste0("(Imputed+NonImputed) log2fc distribution\n n = ", dim(dflog2fc)[1]), xlab = "log2FC")
  garbage <- dev.off()
  
  
  # Relationship between conditions
  # Get the number of biological replicas based on the first condition
  theConditions <- unique(dfmq$GROUP_ORIGINAL)
  theFirstCond <- theConditions[2]
  condFirst <- dfmq[which(dfmq$GROUP_ORIGINAL == theFirstCond),]
  theBiologicalReplicas <- unique(condFirst$SUBJECT_ORIGINAL)
  numberBioReplicas <- length(theBiologicalReplicas)
  
  
  ##############################################################################
  ##############################################################################
  # PLOTS
  
  # boxplot of relative abundances
  cat("Printing out: RELATIVE ABUNDANCE PLOTS\n")
  abundancesName <- gsub(".txt", ".relativeABUNDANCE.pdf", log2fc_file)
  abundancesName <- paste0("plot.",abundancesName)
  abundancesName <- paste0(output_dir,"/",abundancesName)
  
  cat(">> PLOTS: ABUNDANCE PLOTS\n")
  pdf(abundancesName)
    artms_plotAbundanceBoxplots(dfmq)
    artms_plotNumberProteinsAbundance(dfmq)
  garbage <- dev.off()
  
  # Reproducibility plots based on normalized abundance
  cat(">> PLOTS: REPRODUCIBILITY PLOTS\n")
  reproName <- gsub(".txt", ".reproducibilityAbundance.pdf", log2fc_file)
  reproName <- paste0("plot.",reproName)
  reproName <- paste0(output_dir,"/",reproName)
  pdf(reproName)
    artms_plotReproducibilityAbundance(dfmq)
  garbage <- dev.off()
  
  # Conditions
  cat("Printing out: relationship between CONDITIONS\n")
  relaCond <- gsub(".txt", ".relationshipConditions.pdf", log2fc_file)
  relaCond <- paste0("plot.",relaCond)
  relaCond <- paste0(output_dir,"/",relaCond)
  pdf(relaCond)
  plotRelationConditions(dfmq, numberBioReplicas)
  garbage <- dev.off()
  
  # Relationship between log2fc comparisons
  if (length(unique(dflog2fc$Label)) > 1){
    cat("Printing out: log2fc relationships PLOTS\n")
    relaChanges <- gsub(".txt", ".relationshipChanges.pdf", log2fc_file)
    relaChanges <- paste0("plot.",relaChanges)
    relaChanges <- paste0(output_dir,"/",relaChanges)
    pdf(relaChanges)
    plotRatioLog2fc(dflog2fc)
    garbage <- dev.off()
  }else{
    cat("Only one Comparison is available\n")
  }
  
  ##############################################################################
  ##############################################################################
  
  
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ##############################################################################
  ##############################################################################
  # ABUNDANCE DATA, to CREATE FILTERS
  abundance <- artms_loadModelqcBasic(dfmq)
  names(abundance)[grep('Protein', names(abundance))] <- 'Prey'
  names(abundance)[grep('Condition', names(abundance))] <- 'Bait'
  # TECHNICAL REPLICAS: if there are technical replicas, this means that we will find
  # two values for the same protein in the same bioreplica, therefore we need to 
  # aggregate first just in case:
  abundance <- aggregate(Abundance~Prey+Bait+Bioreplica, data = abundance, FUN = mean)
  
  # Let's aggregate to get the sum of the abundance, we will use it later.
  abundance_dcsum <- dcast(abundance, Prey~Bait, value.var = 'Abundance', fun.aggregate = sum, fill = 0 )
  abundance_dcmean <- dcast(abundance, Prey~Bait, value.var = 'Abundance', fun.aggregate = mean, fill = 0 )
  
  
  #########################################################
  # HEATMAPs
  # Use the sum for the heatmap
  dchm_input <- abundance_dcsum
  rownames(dchm_input) <- dchm_input$Prey
  dfhm <- subset(dchm_input, select=-c(Prey))
  aqui <- data.matrix(dfhm)
  
  outHeatMapOverall <- gsub(".txt",".clustering.abundance.all-overview.pdf",log2fc_file)
  outHeatMapOverall <- paste0(output_dir,"/",outHeatMapOverall)
  pheatmap(aqui, filename=outHeatMapOverall, cellwidth=20, main = "Clustered Relative Abundance", cluster_cols = F, fontfamily="Helvetica", labels_row = "", fontsize=6, fontsize_row=8, fontsize_col=8, border_color=NA, fontfamily="Helvetica")
  outHeatMapZoom <- gsub(".txt",".clustering.abundance.all-zoom.pdf",log2fc_file)
  outHeatMapZoom <- paste0(output_dir,"/",outHeatMapZoom)
  pheatmap(aqui, filename=outHeatMapZoom, cellheight = 10, cellwidth=20, main = "Clustered Relative Abundance", cluster_cols = F, fontsize=6, fontsize_row=8, fontsize_col=8, border_color=NA, fontfamily="Helvetica")
  # HEATMAPs
  #########################################################
  
  # Melt again the sum and mean
  abundancelongsum <- melt(abundance_dcsum, id.vars = c('Prey'), value.name = 'Abundance', variable.name = 'Bait' )
  abundancelongmean <- melt(abundance_dcmean, id.vars = c('Prey'), value.name = 'Abundance', variable.name = 'Bait' )
  
  # We dont need the 0 values
  abundancelongsum <- abundancelongsum[!(abundancelongsum$Abundance == 0),]
  abundancelongmean <- abundancelongmean[!(abundancelongmean$Abundance == 0),]
  # Rename and merge:
  names(abundancelongsum)[grep('Abundance', names(abundancelongsum))] <- 'AbSum'
  names(abundancelongmean)[grep('Abundance', names(abundancelongmean))] <- 'AbMean'
  
  abundancelongsummean <- merge(abundancelongsum, abundancelongmean, by = c('Prey', 'Bait'))
  
  # More REPRODUCIBILITY AND SPECIFICY PARATEMERS
  
  # Now the number of bioreplicas based on abundance data
  abundance_dc_length <- dcast(abundance, Prey~Bait, value.var = 'Abundance', fun.aggregate = length, fill = 0 )
  abundancelong_len <- melt(abundance_dc_length, id.vars = c('Prey'), value.name = 'Abundance', variable.name = 'Bait' )
  abundancelong_len <- abundancelong_len[!(abundancelong_len$Abundance == 0),]
  names(abundancelong_len)[grep('Abundance', names(abundancelong_len))] <- 'BioRep'
  
  # IT SHOULD BE A FUNCTION FROM HERE: CALCULATE THE SUM COUNTS From MSrepro
  # Let's dcast and aggregate by condition
  OUTreprod <- dcast(data = abundancelong_len, Prey~Bait, value.var = 'BioRep')
  here <- dim(OUTreprod)[2]
  OUTreprod[is.na(OUTreprod)] <- 0
  # Make a copy to use later
  bioReplicaInfo <- OUTreprod 
  OUTreprod$ReproBioreplicaCount <- rowSums(OUTreprod[,2:here])
  
  reprospec2merge <- subset(OUTreprod, select = c(Prey, ReproBioreplicaCount))
  
  OUTreproCondition <- dcast(data = abundancelong_len, Prey~Bait, value.var = 'BioRep')
  thedim <- dim(OUTreproCondition)[2]
  OUTreproCondition[is.na(OUTreproCondition)] <- 0
  thepreys <- subset(OUTreproCondition, select = c(Prey))
  thevalues <- subset(OUTreproCondition, select = -c(Prey))
  thevalues[thevalues > 0] <- 1
  thevalues$ReproConditionCount <- rowSums(thevalues)
  FinalReproCondition <- cbind(thepreys,thevalues)
  
  reprocondition2merge <- subset(FinalReproCondition, select = c(Prey, ReproConditionCount))
  
  # This version will be printed out below
  OUTreprodFinal <- merge(OUTreprod, reprocondition2merge, by = 'Prey')
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  ####################################################################################
  ####################################################################################
  ## PCA ANALYSIS
  # It requires a simplified version fo modelqc
  # Now let's add the annotation (although I only need the gene name)
  modelqcabundance <- loadModelQCstrict(dfmq, specie)
  cat("PRINCIPAL COMPONENT CHARTS....\n")
  out.pca <- gsub(".txt", "-pca", log2fc_file)
  out.pca <- paste0(output_dir,"/",out.pca)
  getPCAplots(modelqcabundance, out.pca, conditions)
  cat("PCA DONE!\n")
  ####################################################################################
  ####################################################################################
  
  
  modelqc_file_splc <- loadModelQC(dfmq, abundance_dc_length, specie)
  # Now get ready for annotation
  if( grepl("yesptm",isPtm) ){
    names(modelqc_file_splc)[grep('^Protein$', names(modelqc_file_splc))] <- 'Uniprot_PTM'
    # Take the Protein ID, but being very careful about the fluomics labeling
    modelqc_file_splc$Protein <- ifelse(grepl("_H1N1|_H3N2|_H5N1", modelqc_file_splc$Uniprot_PTM), gsub("^(\\S+?_H[1,3,5]N[1,2])_.*", "\\1", modelqc_file_splc$Uniprot_PTM, perl = T) , gsub("^(\\S+?)_.*", "\\1", modelqc_file_splc$Uniprot_PTM, perl = T)) 
    modelqc_file_splc <- artms_annotationUniprot(modelqc_file_splc, 'Protein', specie)
  }else{
    modelqc_file_splc <- artms_annotationUniprot(modelqc_file_splc, 'Protein', specie)
  }
  
  log2fc_file_splc <- loadL2FCWide(dflog2fc, abundance_dc_length, specie)
  # Now get ready for annotation
  if( grepl("yesptm",isPtm) ){
    names(log2fc_file_splc)[grep('^Protein$', names(log2fc_file_splc))] <- 'Uniprot_PTM'
    # Take the Protein ID, but being very careful about the fluomics labeling
    log2fc_file_splc$Protein <- ifelse(grepl("_H1N1|_H3N2|_H5N1", log2fc_file_splc$Uniprot_PTM), gsub("^(\\S+?_H[1,3,5]N[1,2])_.*", "\\1", log2fc_file_splc$Uniprot_PTM, perl = T) , gsub("^(\\S+?)_.*", "\\1", log2fc_file_splc$Uniprot_PTM, perl = T)) 
    log2fc_file_splc <- artms_annotationUniprot(log2fc_file_splc, 'Protein', specie)
  }else{
    log2fc_file_splc <- artms_annotationUniprot(log2fc_file_splc, 'Protein', specie)
  }
  
  log2fc_long <- loadL2FCLong(dflog2fc, specie)
  
  # Now get ready for annotation
  if( grepl("yesptm",isPtm) ){
    names(log2fc_long)[grep('^Protein$', names(log2fc_long))] <- 'Uniprot_PTM'
    # Take the Protein ID, but being very careful about the fluomics labeling
    log2fc_long$Protein <- ifelse(grepl("_H1N1|_H3N2|_H5N1", log2fc_long$Uniprot_PTM), gsub("^(\\S+?_H[1,3,5]N[1,2])_.*", "\\1", log2fc_long$Uniprot_PTM, perl = T) , gsub("^(\\S+?)_.*", "\\1", log2fc_long$Uniprot_PTM, perl = T)) 
    log2fc_long <- artms_annotationUniprot(log2fc_long, 'Protein', specie)
  }else{
    log2fc_long <- artms_annotationUniprot(log2fc_long, 'Protein', specie)
  }
  
  
  # # THE JITTER PLOTS for log2fc values
  # log2fc_long$Specie <- ifelse(grepl("H1N1|H3N2|H5N1", log2fc_long$Protein), "Influenza", "Human")
  # ggplot(log2fc_long[which(log2fc_long$log2FC > -10 & log2fc_long$log2FC < 10),] %>% arrange(Specie), aes(Comparison,log2FC)) +
  #   geom_jitter(aes(colour = Specie), width = 0.5) +
  #   theme_minimal()+
  #   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
  
  
  ################################################################################
  ################################################################################
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # FILTERING BASED ON NUMBER OF BIOLOGICAL REPLICAS
  
  cat("FILTERING BY NUMBER OF BIOLOGICAL REPLICAS\n")
  
  # # This would be the old imputation method
  # imputingMaxMin <- imputingOnTheMaximum(dflog2fc)
  
  # FILTER BASED ON NUMBER OF BIOLOGICAL REPLICAS
  imputedDF <- dflog2fc[c('Protein', 'Label', 'log2FC', 'pvalue', 'adj.pvalue','imputed','iLog2FC','iPvalue')]
  
  # Merge with bioReplicaInf
  cat("Merging with bioReplica Info\n")
  imputedDF <- merge(imputedDF, bioReplicaInfo, by.x = 'Protein', by.y = 'Prey', all.x = T)
  cat("Removing NA\n")
  imputedDF <- imputedDF[!is.na(imputedDF$log2FC),]
  
  # Remove those entries without the minimal number of biological replicas
  # mnbr <- 2 # Minimal Number of Biological Replicas
  cat("Add info about condition more abundant\n")
  imputedDF$CMA <- mapply(selectTheOneLog2fc, imputedDF$iLog2FC, imputedDF$Label)
  
  # this will be the loop for each condition
  theComparisons2check <- unique(imputedDF$Label)
  
  cat("Loop to remove based on number of biological replicas\n")
  for (onlyonecomp in (theComparisons2check)){
    
    ax <- gsub("(.*)(-)(.*)", "\\1", onlyonecomp )
    ay <- gsub("(.*)(-)(.*)", "\\3", onlyonecomp )
    # Remove if does not meet the minimal number of biological replicas in at least one of the conditions, remove
    # BUT CAUTION, it has to be in ONLY the condition
    # Old stuff
    # imputedDF <- imputedDF[-which(((imputedDF[[ax]] < mnbr) & (imputedDF[[ay]] < mnbr)) & (imputedDF$Label == onlyonecomp)),]
    # imputedDF <- imputedDF[-which( (imputedDF[[ax]] < mnbr) & (imputedDF[[ay]] < mnbr) & (imputedDF$Label == onlyonecomp) ),]
    
    # New code: if the condition is not met, i.e., if all the proteins are found in at least X biological replicas, then
    # it would remove the whole thing. 
    if( dim( imputedDF[which( (imputedDF[[ax]] < mnbr) & (imputedDF[[ay]] < mnbr) ),] )[1] > 0){
      imputedDF <- imputedDF[-which( ((imputedDF[[ax]] < mnbr) & (imputedDF[[ay]] < mnbr)) & (imputedDF$Label == onlyonecomp) ),]
    }
    
  }
  cat("Filtering is done!\n")
  
  distributionsFilteredOut <- gsub(".txt",".distributionsFil.pdf",log2fc_file)
  distributionsFilteredOut <- paste0(output_dir,"/",distributionsFilteredOut)
  pdf(distributionsFilteredOut)
  hist(imputedDF$iLog2FC, breaks = 1000, main = paste0("Filtered Log2FC (>2BR)\n n = ",dim(imputedDF)[1]), xlab = "log2fc")
  hist(imputedDF$iPvalue, breaks = 1000, main = paste0("Filtered p-values (>2BR)\n n = ",dim(imputedDF)[1]), xlab = "p-value")
  garbage <- dev.off()
  
  # Stats about imputed values
  yesimputed <- dim(imputedDF[which(imputedDF$imputed=='yes'),])[1]
  nonimputed <- dim(imputedDF[which(imputedDF$imputed=='no'),])[1]
  
  dat <- data.frame(count=c(yesimputed, nonimputed), category = c("Imputed", "Non-Imputed"))
  # Add addition columns, needed for drawing with geom_rect.
  dat$fraction = dat$count / sum(dat$count)
  dat = dat[order(dat$fraction), ]
  dat$ymax = cumsum(dat$fraction)
  dat$ymin = c(0, head(dat$ymax, n=-1))
  
  # ppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
  p1 <- ggplot(dat, aes(fill=category, ymax=ymax, ymin=ymin, xmax=4, xmin=3)) +
    geom_rect() +
    coord_polar(theta="y") +
    xlim(c(0, 4)) +
    labs(title="Proportion of Imputed Intensity values")
  p2 <- ggplot(dat, aes(x=category, y=count, fill=category)) +
    geom_bar(stat="identity") + 
    labs(title="Proportion of Imputed Intensity values")
  
  outImputation <- gsub(".txt",".imputation.pdf",log2fc_file)
  outImputation <- paste0(output_dir,"/",outImputation)
  pdf(outImputation)
  print(p1)
  print(p2)
  garbage <- dev.off()
  # ppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
  
  # End of Imputation
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ################################################################################
  ################################################################################
  
  l2fcol <- dcast(data=imputedDF, Protein~Label, value.var = 'iLog2FC')
  
  ################################################################################
  ################################################################################
  # All the log2fc values
  cat("HEATMAPS with imputed log2fc values\n")
  rownames(l2fcol) <- l2fcol$Protein
  l2fcol <- within(l2fcol, rm(Protein))
  l2fcol[is.na(l2fcol)] <- 0
  
  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # HEATMAPS from log2fc values
  # NA values can be changed to 0
  # Inf values to the maximum value in that comparison
  # -Inf values, changing it to the lowest value
  # Of course, at least 2 comparisons are required to do the clustering
  if(numberConditions > 1){
    l2fcolmatrix <- data.matrix(l2fcol)
    outHeatMapOverallL2fc <- gsub(".txt",".clustering.log2fc.all-overview.pdf",log2fc_file)
    outHeatMapOverallL2fc <- paste0(output_dir,"/",outHeatMapOverallL2fc)
    pheatmap(l2fcolmatrix, filename=outHeatMapOverallL2fc, cellwidth=20, main = "Clustering Log2FC", cluster_cols = F, clustering_method = "average", fontfamily="Helvetica", show_colnames = F, fontsize=6, fontsize_row=3, fontsize_col=10, border_color=NA)
    outHeatMapZoomL2fc <- gsub(".txt",".clustering.log2fc.all-zoom.pdf",log2fc_file)
    outHeatMapZoomL2fc <- paste0(output_dir,"/",outHeatMapZoomL2fc)
    pheatmap(l2fcolmatrix, filename=outHeatMapZoomL2fc, cellheight = 10, cellwidth=20, main = "Clustering Log2FC", cluster_cols = F, fontsize=6, fontsize_row=8, fontsize_col=8, border_color=NA, fontfamily="Helvetica")
    
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # Only significant pvalues
    
    imputedDFsig <- imputedDF[which(imputedDF$iPvalue < 0.05),]
    l2fcolSignificants <- dcast(data=imputedDFsig, Protein~Label, value.var = 'iLog2FC')
    rownames(l2fcolSignificants) <- l2fcolSignificants$Protein
    l2fcolSignificants <- within(l2fcolSignificants, rm(Protein))
    l2fcolSignificants[is.na(l2fcolSignificants)] <- 0
    
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # HEATMAPS from log2fc values
    l2fcolSignificantsmatrix <- data.matrix(l2fcolSignificants)
    outHeatMapOverallL2fc <- gsub(".txt",".clustering.log2fcSign.all-overview.pdf",log2fc_file)
    outHeatMapOverallL2fc <- paste0(output_dir,"/",outHeatMapOverallL2fc)
    pheatmap(l2fcolSignificantsmatrix, filename=outHeatMapOverallL2fc, cellwidth=20, main = "Clustering Log2FC (p-value < 0.05)", cluster_cols = F, fontfamily="Helvetica", labels_row = "", fontsize=6, fontsize_row=8, fontsize_col=8, border_color=NA, fontfamily="Helvetica")
    outHeatMapZoomL2fc <- gsub(".txt",".clustering.log2fcSign.all-zoom.pdf",log2fc_file)
    outHeatMapZoomL2fc <- paste0(output_dir,"/",outHeatMapZoomL2fc)
    pheatmap(l2fcolSignificantsmatrix, filename=outHeatMapZoomL2fc, cellheight = 10, cellwidth=20, main = "Clustering Log2FC (p-value < 0.05)", cluster_cols = F, fontsize=6, fontsize_row=8, fontsize_col=8, border_color=NA, fontfamily="Helvetica")
    cat("...and heatmaps are done\n")
  }
  ################################################################################
  
  
  # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  # ENRICHMENT OF MOST ABUNDANT PROTEINS (from IMPUTED LOG2FC values)
  # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  # Let's select first significance based on pvalue, by using the iPvalue we are already including the imputed pvalues...
  l2fcolcopy <- dcast(data=imputedDF[which(imputedDF$iPvalue < 0.05),], Protein~Label, value.var = 'iLog2FC')
  
  if(dim(l2fcolcopy)[1] > 0){
    # Let's melt now for enrichment analysis
    l2fcol4enrichment <- melt(data=l2fcolcopy, id.vars = c('Protein'))
    names(l2fcol4enrichment)[grep('variable', names(l2fcol4enrichment))] <- 'Comparisons'
    l2fcol4enrichment <- artms_annotationUniprot(l2fcol4enrichment, 'Protein', specie)
    l2fcol4enrichment <- l2fcol4enrichment[!is.na(l2fcol4enrichment$value),]
  }
  
  if(enrich == "yesenrich" & dim(l2fcolcopy)[1] > 0 ){
    
    # Now it is chosen by the user
    # threshold <- 1
    
    if( grepl("yesptm",isPtm) ){
      l2fcol4enrichment <- within(l2fcol4enrichment, rm(Gene,Entry.name,Protein.names))
      # Remove parties for enrichment
      l2fcol4enrichment <- l2fcol4enrichment[grep(",",l2fcol4enrichment$Protein, invert=T),]
      # Select the Uniprot ID, but keep in mind that some of them might have many _ph54_ph446
      # before
      # l2fcol4enrichment$Protein <- gsub("(.*)(_.*)", "\\1", l2fcol4enrichment$Protein)
      l2fcol4enrichment$Protein <- gsub("^(\\S+?)_.*", "\\1", l2fcol4enrichment$Protein, perl = T)
      l2fcol4enrichment <- unique(l2fcol4enrichment)
      l2fcol4enrichment <- proteinid2gene(l2fcol4enrichment, specie, 'Protein')
    }
    
    # ALL SIGNIFICANT CHANGES log2fc
    # GPROFILER
    cat("GPROFILER Enrichment (ALL significant Changes)\n")
    
    filallsig_log2fc_long <- l2fcol4enrichment[which(abs(l2fcol4enrichment$value) >= threshold), ]
    
    if(dim(filallsig_log2fc_long)[1] > 0){
      
      out.mac.allsig <- gsub(".txt","-enrich-MAC-allsignificants.txt",log2fc_file)
      out.mac.allsig <- paste0(output_dir,"/",out.mac.allsig)
      
      mac.allsig <- enrichImputedLog2fc(filallsig_log2fc_long, out.mac.allsig, specie, listOfGenes)
      
      if( dim(mac.allsig)[1] > 0 ){
        write.table(mac.allsig, out.mac.allsig, quote = F, sep = "\t", row.names = F, col.names = T)
      }
      
      cat("Corum Complexes Enrichment (allsignificants MACs)\n")
      
      # CORUM
      allsigComplexEnriched <- enrichForComplexes(filallsig_log2fc_long, backgroundNumber)
      
      if(dim(allsigComplexEnriched)[1] > 0){
        out.mac.allsig.corum <- gsub(".txt","-enrich-MAC-allsignificants-corum.txt",log2fc_file)
        out.mac.allsig.corum <- paste0(output_dir,"/",out.mac.allsig.corum)
        write.table(allsigComplexEnriched, out.mac.allsig.corum, quote = F, sep = "\t", row.names = F, col.names = T)
      }
      
      # And the heatmap
      if(dim(allsigComplexEnriched)[1] > 2){
        cat("Plotting all significants corum complexes\n")
        out.mac.allsig.corum.pdf <- gsub(".txt","-enrich-MAC-allsignificants-corum.pdf",log2fc_file)
        out.mac.allsig.corum.pdf <- paste0(output_dir,"/",out.mac.allsig.corum.pdf)
        # out.mac.allsig.corum.pdf <- 'whatever.corum.allsigitive.pdf'
        plot.corum(allsigComplexEnriched, out.mac.allsig.corum.pdf, "MAC ALL SIGNIFICANT Protein Complex Enrichment")
      }else{
        
        cat("Not enough negative corum complexes to plot\n")
        
      }  
    }else{
      stop("\n\t---- NOTHING is significant! Check what's going on\n\n")
      mac.allsig <- NULL
      allsigComplexEnriched <- NULL
    }
    
    
    # POSITIVE log2fc
    # GPROFILER
    cat("GPROFILER Enrichment (Positive MACs)\n")
    
    filpos_log2fc_long <- l2fcol4enrichment[which(l2fcol4enrichment$value >= threshold), ]
    
    if(dim(filpos_log2fc_long)[1] > 0){
      
      out.mac.pos <- gsub(".txt","-enrich-MAC-positives.txt",log2fc_file)
      out.mac.pos <- paste0(output_dir,"/",out.mac.pos)
      mac.pos <- enrichImputedLog2fc(filpos_log2fc_long, out.mac.pos, specie, listOfGenes)
      
      if(dim(mac.pos)[1] > 0){
        write.table(mac.pos, out.mac.pos, quote = F, sep = "\t", row.names = F, col.names = T)  
      }
      
      cat("Corum Complexes Enrichment (Positives MACs)\n")
      
      # CORUM
      positiveComplexEnriched <- enrichForComplexes(filpos_log2fc_long, backgroundNumber)
      
      if(dim(positiveComplexEnriched)[1] > 0){
        out.mac.pos.corum <- gsub(".txt","-enrich-MAC-positives-corum.txt",log2fc_file)
        out.mac.pos.corum <- paste0(output_dir,"/",out.mac.pos.corum)
        write.table(positiveComplexEnriched, out.mac.pos.corum, quote = F, sep = "\t", row.names = F, col.names = T)
      }
      
      # And the heatmap
      if(dim(positiveComplexEnriched)[1] > 2){
        cat("Plotting positive corum complexes\n")
        out.mac.pos.corum.pdf <- gsub(".txt","-enrich-MAC-positives-corum.pdf",log2fc_file)
        out.mac.pos.corum.pdf <- paste0(output_dir,"/",out.mac.pos.corum.pdf)
        # out.mac.pos.corum.pdf <- 'whatever.corum.positive.pdf'
        plot.corum(positiveComplexEnriched, out.mac.pos.corum.pdf, "MAC+ Protein Complex Enrichment")
      }else{
        
        cat("Not enough negative corum complexes to plot\n")
        
      }    
    }else{
      cat("\t\t ------ Nothing is significant in the Positive site of things")
      mac.pos <- NULL
      positiveComplexEnriched <- NULL
    }
    
    
    # NEGATIVE log2fc
    cat("GPROFILER Enrichment (Negative MACs)\n")
    
    filneg_log2fc_long <- l2fcol4enrichment[which(l2fcol4enrichment$value <= -threshold), ]
    
    if(dim(filneg_log2fc_long)[1] > 0){
      
      out.mac.neg <- gsub(".txt","-enrich-MAC-negatives.txt",log2fc_file)
      out.mac.neg <- paste0(output_dir,"/",out.mac.neg)  
      mac.neg <- enrichImputedLog2fc(filneg_log2fc_long, out.mac.neg, specie, listOfGenes)
      
      if( dim(mac.neg)[1] > 0 ){
        write.table(mac.neg, out.mac.neg, quote = F, sep = "\t", row.names = F, col.names = T)
      }
      
      cat("Corum Complexes Enrichment (Negative MACs)\n")
      
      negativesComplexEnriched <- enrichForComplexes(filneg_log2fc_long, backgroundNumber)
      
      if(dim(negativesComplexEnriched)[1] > 0){
        out.mac.neg.corum <- gsub(".txt","-enrich-MAC-negatives-corum.txt",log2fc_file)
        out.mac.neg.corum <- paste0(output_dir,"/",out.mac.neg.corum)
        write.table(negativesComplexEnriched, out.mac.neg.corum, quote = F, sep = "\t", row.names = F, col.names = T)  
      }
      
      # And the heatmap
      if(dim(negativesComplexEnriched)[1] > 2){
        cat("Plotting negative corum complexes\n")
        out.mac.neg.corum.pdf <- gsub(".txt","-enrich-MAC-negatives-corum.pdf",log2fc_file)
        out.mac.neg.corum.pdf <- paste0(output_dir,"/",out.mac.neg.corum.pdf)
        plot.corum(negativesComplexEnriched, out.mac.neg.corum.pdf, "MAC- Protein Complex Enrichment")
      }else{
        
        cat("Not enough negative corum complexes to plot\n")
        
      }    
    }else{
      cat("\t\t ------ Nothing is significant in the Positive site of things")
      mac.neg <- NULL
      negativesComplexEnriched <- NULL
    }
  }else{
    cat("NO ENRICHMENT of log2fc values chosen\n")
  }
  # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  superunified <- merge(abundancelongsummean, abundancelong_len, by = c('Bait', 'Prey'), all = T)
  superunified <- merge(superunified, reprocondition2merge, by = 'Prey', all = T)
  superunified <- merge(superunified, reprospec2merge, by = 'Prey', all = T)
  
  if( grepl("yesptm",isPtm) ){
    names(superunified)[grep('^Prey$', names(superunified))] <- 'Uniprot_PTM'
    # Take the Protein ID, but being very careful about the fluomics labeling
    superunified$Prey <- ifelse(grepl("_H1N1|_H3N2|_H5N1", superunified$Uniprot_PTM), gsub("^(\\S+?_H[1,3,5]N[1,2])_.*", "\\1", superunified$Uniprot_PTM, perl = T) , gsub("^(\\S+?)_.*", "\\1", superunified$Uniprot_PTM, perl = T)) 
  }
  
  superunified <- artms_annotationUniprot(superunified, 'Prey', specie)
  
  # Rename (before it was just a lazy way to use another code)
  names(superunified)[grep('Bait', names(superunified))] <- 'Condition'
  
  # ANNOTATE SPECIE
  cat("Annotating the specie before jitter plot\n")
  
  annotateSpecie <- function(df, pathogen, specie){
    if(pathogen == "nopathogen"){
      # Influenza is treated differently
      df$Specie <- ifelse(grepl("_H1N1|_H3N2|_H5N1", df$Protein), "Influenza", specie)  
    }else{
      # Pathogens
      df$Specie <- ifelse(df$Protein %in% pathogen.ids$Entry, pathogen, specie)
    }
    return(df)
  }
  
  superunified <- annotateSpecie(superunified, pathogen, specie)
  imputedDF <- annotateSpecie(imputedDF, pathogen, specie)
  
  
  
  #########################################################
  #########################################################
  # THE JITTER PLOTS
  # TO DO: REMOVE PROTEINS NOT FOUND IN MORE THAN 1 BIOLOGICAL REPLICA
  cat("Time for the Jitter plot\n")
  if(isFluomics == "yesflu"){
    # Filter by number of biological replicas > 1
    superunifiedfiltered <- superunified[which(superunified$BioRep > 1),]
    # Filter less than 
    superunifiedfiltered <- superunifiedfiltered[which(superunifiedfiltered$AbMean < 30 & superunifiedfiltered$AbMean > 10 ), ]
    # Removing carry overs
    superunifiedfiltered <- superunifiedfiltered[!(grepl("H1N1|H3N2|H5N1", superunifiedfiltered$Protein) & grepl("MOCK",superunifiedfiltered$Condition)),]
    
    cat("\n\t>>> Printing out: Jittered PLOTS")
    abuJittered <- gsub(".txt", ".abundanceGrouped.pdf", log2fc_file)
    abuJittered <- paste0("plot.",abuJittered)
    abuJittered <- paste0(output_dir,"/",abuJittered)
    suppressMessages(library(dplyr))
    # j <- ggplot(superunifiedfiltered %>% arrange(Specie), aes(Condition,AbMean))
    # j <- j + geom_jitter(aes(colour = Specie), width = 0.3)
    if(specie == "human"){
      j <- ggplot(superunifiedfiltered %>% arrange(desc(Specie)), aes(x = Condition, y = AbMean, colour = Specie)) #superunifiedfiltered %>% arrange(Specie)
      j <- j + geom_jitter(width = 0.3)
      j <- j + scale_colour_manual(values = c("red", "lightblue"))
    } else if(specie == "mouse"){
      j <- ggplot(superunifiedfiltered %>% arrange(Specie), aes(Condition, AbMean))
      j <- j + geom_jitter(aes(colour = Specie), width = 0.3)
      j <- j + scale_colour_manual(values = c("azure3","red"))
    }
    j <- j + theme_minimal()
    j <- j + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
    
    # Do you want to add labels?
    # if( grepl("yesptm", isPtm) ){
    #    j <- j + geom_text(aes(label=ifelse(Specie=="Influenza",as.character(Uniprot_PTM),'')),hjust=0,vjust=0,size=2)
    #    j <- j + geom_text(aes(label=ifelse(Specie=="TB",as.character(Uniprot_PTM),'')),hjust=0,vjust=0,size=2)
    # }else{
    #    j <- j + geom_text(aes(label=ifelse(Specie=="Influenza",as.character(Protein),'')),hjust=0,vjust=0,size=2)
    #    j <- j + geom_text(aes(label=ifelse(Specie=="TB",as.character(Protein),'')),hjust=0,vjust=0,size=2)
    # }
    # j <- j + geom_text_repel(aes(label=ifelse(Specie=="Influenza",as.character(Protein),'')), size=2)
    
    pdf(abuJittered)
    print(j)
    garbage <- dev.off()
    cat("\t--->done\n\n")
  }
  #################################################################################
  #################################################################################
  
  # Wide version of imputed values
  cat("Wide format for imputed values (log2fc only)\n")
  dcImputed <- dcast(data = imputedDF, Protein~Label, value.var = "iLog2FC")
  
  
  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # OUTPUT FILES, txt and xlsx, using the "results" file from MSstats as the template name
  outmodeqcLong <- gsub(".txt","-longAbundance.txt", log2fc_file)
  outmodeqcLong <- paste0(output_dir,"/",outmodeqcLong)
  write.table(superunified, outmodeqcLong, quote = F, sep = "\t", row.names = F, col.names = T)
  
  outmodelqc <- gsub(".txt","-wideAbundance.txt", log2fc_file)
  outmodelqc <- paste0(output_dir,"/",outmodelqc)
  write.table(modelqc_file_splc, outmodelqc, quote = F, sep = "\t", row.names = F, col.names = T)
  
  outlog2fc <- gsub(".txt","-wideL2fc.txt", log2fc_file)
  outlog2fc <- paste0(output_dir,"/",outlog2fc)
  write.table(log2fc_file_splc, outlog2fc, quote = F, sep = "\t", row.names = F, col.names = T)
  
  outwideimputed <- gsub(".txt","-wideImputedL2fc.txt", log2fc_file)
  outwideimputed <- paste0(output_dir,"/",outwideimputed)
  write.table(dcImputed, outwideimputed, quote = F, sep = "\t", row.names = F, col.names = T)
  
  # iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
  # IMPUTED extended: spliting the site information for other purposes
  # Testing to generate the input file for pedro beltrao
  # library(tidyr)
  # library(dplyr)
  
  if( isPtm == "yesptmph" ){
    imputedDFext <- imputedDF
    names(imputedDFext)[grep('^Protein$', names(imputedDFext))] <- 'Uniprot_PTM'
    # Take the Protein ID, but being very careful about the fluomics labeling
    imputedDFext$Protein <- ifelse(
      grepl("_H1N1|_H3N2|_H5N1", imputedDFext$Uniprot_PTM), 
      gsub("^(\\S+?_H[1,3,5]N[1,2])_.*", "\\1", imputedDFext$Uniprot_PTM, perl = T) , 
      gsub("^(\\S+?)_.*", "\\1", imputedDFext$Uniprot_PTM, perl = T)
    ) 
    
    # Extract sites from Uniprot_PTM
    imputedDFext$PTMsite <- gsub("^(\\S+?)(_ph.*)", "\\2", imputedDFext$Uniprot_PTM, perl = T)
    imputedDFext$PTMsite <- gsub("^(_ph)","", imputedDFext$PTMsite)
    imputedDFext$PTMsite <- gsub("_ph", ",", imputedDFext$PTMsite)
    # And create independent columns for each of them
    imputedDFext <- imputedDFext %>% mutate(PTMsite = strsplit(PTMsite, ",")) %>% unnest(PTMsite)
    imputedDFext <- artms_annotationUniprot(imputedDFext, 'Protein', specie)
    names(imputedDFext)[grep("^Label$", names(imputedDFext))] <- 'Comparison'
    
    # to delete
    # imputedDFext$Specie <- ifelse(grepl("_H1N1|_H3N2|_H5N1", imputedDFext$Protein), "Influenza", specie)  
    # imputedDFext$Specie <- ifelse(imputedDFext$Protein %in% pathogen.ids$Entry, pathogen, specie)
    imputedDFext <- annotateSpecie(imputedDFext, pathogen, specie)
    
    outlog2fcImputext <- gsub(".txt","-imputedL2fcExtended.txt", log2fc_file)
    outlog2fcImputext <- paste0(output_dir,"/",outlog2fcImputext)
    write.table(imputedDFext, outlog2fcImputext, quote = F, sep = "\t", row.names = F, col.names = T)
  } else if(isPtm == "yesptmsites") {
    imputedDFext <- imputedDF
    #1. Change the Protein name
    names(imputedDFext)[grep('^Protein$', names(imputedDFext))] <- 'Uniprot_PTM'
    
    # 2. Make a copy of Uniprot_PTM to operate on it
    imputedDFext$PTMone <- imputedDFext$Uniprot_PTM
    
    # 3. Create independent columns for each of them
    imputedDFext <- imputedDFext %>% mutate(PTMone = strsplit(PTMone, ",")) %>% unnest(PTMone)
    
    # 4. And take the labels:
    imputedDFext$Protein <- ifelse(grepl("_H1N1|_H3N2|_H5N1", imputedDFext$PTMone), 
                                   gsub("^(\\S+?_H[1,3,5]N[1,2])_.*", "\\1", imputedDFext$PTMone, perl = T) , 
                                   gsub("^(\\S+?)_.*", "\\1", imputedDFext$PTMone, perl = T)) 
    imputedDFext$PTMsite <- gsub("(\\S+)(_[S,T,Y,K])(\\d+)","\\3",imputedDFext$PTMone)
    
    imputedDFext <- artms_annotationUniprot(imputedDFext, 'Protein', specie)
    names(imputedDFext)[grep("^Label$", names(imputedDFext))] <- 'Comparison'
    
    # imputedDFext$Specie <- ifelse(grepl("_H1N1|_H3N2|_H5N1", imputedDFext$Protein), "Influenza", specie)  
    # imputedDFext$Specie <- ifelse(imputedDFext$Protein %in% pathogen.ids$Entry, pathogen, specie)  
    imputedDFext <- annotateSpecie(imputedDFext, pathogen, specie)
    
    outlog2fcImputext <- gsub(".txt","-imputedL2fcExtended.txt", log2fc_file)
    outlog2fcImputext <- paste0(output_dir,"/",outlog2fcImputext)
    write.table(imputedDFext, outlog2fcImputext, quote = F, sep = "\t", row.names = F, col.names = T)
  }
  
  # IMPUTED non extended version
  # Let's annotate Imputed table
  
  if( grepl("yesptm",isPtm) ){
    names(imputedDF)[grep('Protein', names(imputedDF))] <- 'Uniprot_PTM'
    imputedDF$UniprotID <- imputedDF$Uniprot_PTM
    # The virus labeling has to be taken into account when getting the uniprot id:
    imputedDF$UniprotID <- ifelse(grepl("_H1N1|_H3N2|_H5N1", imputedDF$UniprotID), gsub("^(\\S+?_H[1,3,5]N[1,2])_.*", "\\1", imputedDF$UniprotID, perl = T) , gsub("^(\\S+?)_.*", "\\1", imputedDF$UniprotID, perl = T)) 
    imputedDF <- artms_annotationUniprot(imputedDF, 'UniprotID', specie)
    names(imputedDF)[grep("Label", names(imputedDF))] <- 'Comparison'
    
    # to delete
    # imputedDF$Specie <- ifelse(imputedDF$Protein %in% pathogen.ids$Entry, pathogen, specie)
    imputedDF <- annotateSpecie(imputedDF, pathogen, specie)
    
    # Wide version of imputedDF
    imputedDF_wide_log2fc <- dcast(data = imputedDF, Gene+Protein+Uniprot_PTM~Comparison, value.var = 'iLog2FC', fill = 0)
    imputedDF_wide_pvalue <- dcast(data = imputedDF, Gene+Protein+Uniprot_PTM~Comparison, value.var = 'iPvalue', fill = 0)
    
  }else if(isPtm == "noptmsites"){
    imputedDF <- artms_annotationUniprot(imputedDF, 'Protein', specie)
    names(imputedDF)[grep("Label", names(imputedDF))] <- 'Comparison'
    
    # imputedDF$Specie <- ifelse(imputedDF$Protein %in% pathogen.ids$Entry, pathogen, specie)
    imputedDF <- annotateSpecie(imputedDF, pathogen, specie)
    
    # Wide version of imputedDF
    imputedDF_wide_log2fc <- dcast(data = imputedDF, Gene+Protein~Comparison, value.var = 'iLog2FC', fill = 0)
    imputedDF_wide_pvalue <- dcast(data = imputedDF, Gene+Protein~Comparison, value.var = 'iPvalue', fill = 0)
  }else{
    stop("you should not see this message. Just go to the source code if you did\n")
  }
  
  
  # boxplot of relative abundances
  cat("Printing final number in imputed Conditions\n")
  numimputedfinal <- gsub(".txt", ".FNIC.pdf", log2fc_file)
  numimputedfinal <- paste0("plot.",numimputedfinal)
  numimputedfinal <- paste0(output_dir,"/",numimputedfinal)
  
  pdf(numimputedfinal)
  plotNumberProteinsImputedLog2fc(imputedDF)
  garbage <- dev.off()
  
  
  # And print it out
  outlog2fcImpute <- gsub(".txt","-imputedL2fc.txt", log2fc_file)
  outlog2fcImpute <- paste0(output_dir,"/",outlog2fcImpute)
  write.table(imputedDF, outlog2fcImpute, quote = F, sep = "\t", row.names = F, col.names = T)
  
  # PCA AND CLUSTERING ANALYSIS
  if(isPtm == "noptmsites"){
    
    cat("The Clustering analysis begins... \n")
    # GET THE LIST OF SIGNIFICANTS FOR THE EXPERIMENT(S)
    list_of_significants <- unique(imputedDF$Protein[which( abs(imputedDF$iLog2FC > 1) & imputedDF$iPvalue < 0.05 )])
    
    # AND APPLY THE FILTER
    data.select <- imputedDF[which(imputedDF$Protein %in% list_of_significants),]
    
    # Two options....
    
    # GENE BASED -> heatmaps
    hasdc <- dcast(data = data.select[which(data.select$imputed == "no"),], Gene+Protein~Comparison, value.var = "iLog2FC", fun.aggregate = median, fill = 0)
    
    # EXPERIMENT BASED -> PCA
    hasdcexp <- dcast(data = data.select[which(data.select$imputed == "no"),], Comparison~Gene+Protein, value.var = "iLog2FC", fun.aggregate = median, fill = 0)
    
    # ----------------------------------------------------------------------
    # CLUSTERING ANALYSIS
    
    # GENE BASED
    rownames(hasdc) <- paste0(hasdc$Gene,"_",hasdc$Protein)
    vamos <- within(hasdc, rm(Gene,Protein))
    venga <- as.matrix(vamos)
    
    # EXPERIMENT BASED
    rownames(hasdcexp) <- hasdcexp$Comparison
    vamosexp <- within(hasdcexp, rm(Comparison))
    vengaexp <- as.matrix(vamosexp)
    
    
    # PCA AND CORRELATION ANALYSIS
    suppressMessages(library("FactoMineR"))
    suppressMessages(library("devtools"))
    suppressMessages(library("factoextra"))
    suppressMessages(library("corrplot"))
    suppressMessages(library("PerformanceAnalytics"))
    
    
    # Correlation matrix
    df.cor.matrix <- round(cor(venga, use = "pairwise.complete.obs"), 2)
    
    file_corr_l2fc <- gsub(".txt",".log2fc-corr.pdf",log2fc_file)
    file_corr_l2fc <- paste0(output_dir,"/",file_corr_l2fc)
    pdf(file_corr_l2fc, width = 12, height = 9)
    corrplot(df.cor.matrix,
             type = "upper",
             tl.pos = "td",
             method = "circle",
             tl.cex = 0.9,
             tl.col = 'black',
             tl.srt=45,
             # order = "hclust",
             diag = T)
    chart.Correlation(venga, histogram=TRUE, pch=25, main = "Correlation between Comparisons")
    garbage <- dev.off()
    
    # BASED ON GROUPS
    pca.hasdcexp <- PCA(hasdcexp[,-c(1)], scale.unit = FALSE, ncp = 4, graph=FALSE)
    
    pca_all <- fviz_pca_ind(pca.hasdcexp,
                            labelsize = 3,
                            repel = TRUE,
                            habillage = as.factor(hasdcexp$Comparison),
                            addEllipses=F,
                            ellipse.level=0.95)
    
    file_pca_l2fc <- gsub(".txt",".log2fc-pca.pdf",log2fc_file)
    file_pca_l2fc <- paste0(output_dir,"/",file_pca_l2fc)
    
    pdf(file_pca_l2fc, width = 9, height = 7)
    print(pca_all)
    garbage <- dev.off()
    
    # Determine the OPTIMAL NUMBER OF CLUSTERS:
    suppressMessages(library("NbClust"))
    suppressMessages(library("cluster"))
    suppressMessages(library("ComplexHeatmap"))
    suppressMessages(library("circlize"))
    
    # Elbow method
    e1 <- fviz_nbclust(venga, kmeans, method = "wss") +
      geom_vline(xintercept = 4, linetype = 2)+
      labs(subtitle = "kmeans Elbow method")
    e2 <- fviz_nbclust(venga, cluster::pam, method = "wss") +
      geom_vline(xintercept = 4, linetype = 2)+
      labs(subtitle = "PAM Elbow method")
    
    # Silhouette method
    k1 <- fviz_nbclust(venga, kmeans, method = "silhouette")+
      labs(subtitle = "kmeans Silhouette method")
    k2 <- fviz_nbclust(venga, cluster::pam, method = "silhouette")+
      labs(subtitle = "pam Silhouette method")
    
    
    # Create a dendrogram
    # library(factoextra)
    res.dist <- get_dist(vamosexp, stand = TRUE, method = "minkowski")
    hc <- hclust(res.dist)
    file_dendro_l2fc <- gsub(".txt",".log2fc-dendro.pdf",log2fc_file)
    file_dendro_l2fc <- paste0(output_dir,"/",file_dendro_l2fc)
    pdf(file_dendro_l2fc, width = 9, height = 7)
    plot(hc)
    garbage <- dev.off()
    
    # COMPLEXHEATMAP Heatmap with a specified number of optimal clusters
    n = 10
    pam.res <- pam(vamos, k=n)
    
    cp1 <- fviz_cluster(pam.res)
    cp2 <- fviz_silhouette(silhouette(pam.res))
    
    file_clusterplots_l2fc <- gsub(".txt",".log2fc-clusters.pdf",log2fc_file)
    file_clusterplots_l2fc <- paste0(output_dir,"/",file_clusterplots_l2fc)
    pdf(file_clusterplots_l2fc, width = 9, height = 7)
    print(e1)
    print(e2)
    print(k1)
    print(k2)
    print(cp1)
    print(cp2)
    garbage <- dev.off()
    
    hmap <- Heatmap(vamos,
                    name=paste0("Clusters ","(n = ",n,")"),
                    col = circlize::colorRamp2(c(-3, 0, 3), c("firebrick1", "black", "olivedrab1")),
                    heatmap_legend_param=list(color_bar="continuous", legend_direction="horizontal", legend_width=unit(5,"cm"), title_position="topcenter", title_gp=gpar(fontsize=15, fontface="bold")),
                    split=paste0("", pam.res$clustering),
                    row_title="Genes",
                    row_title_side="left",
                    row_title_gp=gpar(fontsize=15, fontface="bold"),
                    show_row_names=FALSE,
                    column_title="Relative Quantifications",
                    column_title_side="top",
                    column_title_gp=gpar(fontsize=10, fontface="bold"),
                    column_title_rot=0,
                    show_column_names=TRUE,
                    cluster_columns = F,
                    clustering_distance_columns=function(x) as.dist(1-cor(t(x))),
                    clustering_method_columns="ward.D2",
                    clustering_distance_rows="euclidean",
                    clustering_method_rows="ward.D2",
                    row_dend_width=unit(30,"mm"),
                    column_dend_height=unit(30,"mm"),
                    # top_annotation=colAnn,
                    top_annotation_height=unit(1.75,"cm"),
                    # bottom_annotation=sampleBoxplot,
                    bottom_annotation_height=unit(4, "cm"),
                    column_names_gp = gpar(fontsize = 10))
    
    file_clusterheat_l2fc <- gsub(".txt",".log2fc-clusterheatmap.pdf",log2fc_file)
    file_clusterheat_l2fc <- paste0(output_dir,"/",file_clusterheat_l2fc)
    pdf(file_clusterheat_l2fc, width = 12, height = 10)
    draw(hmap, heatmap_legend_side="top", annotation_legend_side="right")
    garbage <- dev.off()
    
    cl_number <- pam.res$clustering
    dfclusters <- as.data.frame(cl_number)
    dfclusters$ids <- row.names(dfclusters)
    dfclusters$Gene <- gsub("(.*)(_)(.*)","\\1",dfclusters$ids)
    dfclusters$Protein <- gsub("(.*)(_)(.*)","\\3",dfclusters$ids)
    
    # Making sure we have unique genes in each comparison (the PTM might bring redundancy)
    pretmp <- dfclusters[c('Gene', 'cl_number')]
    pretmp <- unique(pretmp)
    
    tmp = split(pretmp$Gene, pretmp$cl_number, drop=T)
    
    if(specie == "human"){
      enrichgenes <- enrichProfiler(tmp, categorySource = c('GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC', 'CORUM', 'HPA','OMIM'), specie = 'hsapiens', listOfGenes) # 'HP'
    }else if(specie == "mouse"){
      enrichgenes <- enrichProfiler(tmp, categorySource = c('GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC', 'CORUM'), specie = 'mmusculus', listOfGenes)
    }else{
      stop("\n\n\ntOhhh no, this specie",specie," is not supported in the enrichment!!\n\n\n")
    }
    
    file_clusterheatenrich_l2fc <- gsub(".txt",".log2fc-clusterheatmap-enriched.txt",log2fc_file)
    file_clusterheatenrich_l2fc <- paste0(output_dir,"/",file_clusterheatenrich_l2fc)
    write.table(enrichgenes, file_clusterheatenrich_l2fc, col.names = T, row.names = F, sep = "\t", quote = F)
    
    file_clusterheatdata_l2fc <- gsub(".txt",".log2fc-clusterheatmap.txt",log2fc_file)
    file_clusterheatdata_l2fc <- paste0(output_dir,"/",file_clusterheatdata_l2fc)
    write.table(dfclusters, file_clusterheatdata_l2fc, col.names = T, row.names = F, sep = "\t", quote = F)
  }
  
  
  
  # iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
  
  # Move Gene name to the left:
  # Merge with bioReplicaInf
  cat("Merging with bioReplica Info to long_log2fc\n")
  # This depends on the type of dataset
  if( grepl("yesptm",isPtm) ){
    log2fc_long <- merge(log2fc_long, bioReplicaInfo, by.x = 'Uniprot_PTM', by.y = 'Prey', all.x = T) 
  }else{
    log2fc_long <- merge(log2fc_long, bioReplicaInfo, by.x = 'Protein', by.y = 'Prey', all.x = T)
  }
  
  cat("Removing NA\n")
  log2fc_long <- log2fc_long[!is.na(log2fc_long$log2FC),]
  
  # Merging lo
  outlog2fclong <- gsub(".txt","-longL2fc.txt", log2fc_file)
  outlog2fclong <- paste0(output_dir,"/",outlog2fclong)
  write.table(log2fc_long, outlog2fclong, quote = F, sep = "\t", row.names = F, col.names = T)
  
  # outUniqueProteinsCondition <- gsub(".txt","-uniquePerCondition.txt", log2fc_file)
  # outUniqueProteinsCondition <- paste0(output_dir,"/",outUniqueProteinsCondition)
  # write.table(aocabmerge, outUniqueProteinsCondition, quote = F, sep = "\t", row.names = F, col.names = T)
  
  outexcel <- gsub(".txt","-summary.xlsx",log2fc_file)
  outexcel <- paste0(output_dir,"/",outexcel)
  
  if(enrich == "yesenrich"){
    # But now check whether is a PTM case:
    if( grepl("yesptm",isPtm) ){
      list_of_datasets <- list(
        # "AbundanceLong" = superunified,
        # "AbundanceWide" = modelqc_file_splc,
        # "log2fcWide" = log2fc_file_splc,
        # "log2fcLong" = log2fc_long,
        "log2fcImputed" = imputedDF,
        "log2fcImpExt" = imputedDFext,
        "wide_iLog2fc" = imputedDF_wide_log2fc,
        "wide_iPvalue" = imputedDF_wide_pvalue,
        "enrichALL" = mac.allsig,
        "enrichMACpos" = mac.pos,
        "enrichMACneg" = mac.neg,
        "enMACallCorum" = allsigComplexEnriched,
        "enMACposCorum" = positiveComplexEnriched,
        "enMACnegCorum" = negativesComplexEnriched)
    }else if(isPtm == "noptmsites"){
      list_of_datasets <- list(
        # "AbundanceLong" = superunified,
        # "AbundanceWide" = modelqc_file_splc,
        # "log2fcWide" = log2fc_file_splc, 
        # "log2fcLong" = log2fc_long, 
        "log2fcImputed" = imputedDF,
        "wide_iLog2fc" = imputedDF_wide_log2fc,
        "wide_iPvalue" = imputedDF_wide_pvalue,
        "enrichALL" = mac.allsig,
        "enrich-MACpos" = mac.pos,
        "enrich-MACneg" = mac.neg,
        "enMACallCorum" = allsigComplexEnriched,
        "enMACposCorum" = positiveComplexEnriched,
        "enMACnegCorum" = negativesComplexEnriched)
      
    }else{
      stop("Oh no!! This will fail if you are using UB!!\n")
    }
  }else if(enrich == "noenrich"){
    cat("\t\t-----+ You chose not to enrich\n")
    if( grepl("yesptm",isPtm) ) {
      list_of_datasets <- list(
        # "AbundanceLong" = superunified,
        # "AbundanceWide" = modelqc_file_splc,
        # "log2fcWide" = log2fc_file_splc,
        # "log2fcLong" = log2fc_long,
        "log2fcImputed" = imputedDF,
        "log2fcImpExt" = imputedDFext,
        "wide_iLog2fc" = imputedDF_wide_log2fc,
        "wide_iPvalue" = imputedDF_wide_pvalue)
    }else if(isPtm == "noptmsites"){
      list_of_datasets <- list(
        # "AbundanceLong" = superunified,
        # "AbundanceWide" = modelqc_file_splc,
        # "log2fcWide" = log2fc_file_splc, 
        # "log2fcLong" = log2fc_long, 
        "log2fcImputed" = imputedDF,
        "wide_iLog2fc" = imputedDF_wide_log2fc,
        "wide_iPvalue" = imputedDF_wide_pvalue
      )
    }else{
      stop("Oh no!! This will fail if you are using UB!!\n")
    }
  }else{
    stop("\n\nYOU SHOULD NEVER SEE THIS MESSAGE. IF you do, dude, check the source code urgently\n\n")
    # The script should have crashed by this point. If it gets up to here... it would be very weird
  }
  
  # Defining style for the header
  hs <- createStyle(fontName = "Arial", fontColour = "white", fgFill = "#000000",
                    textDecoration = "Bold", border = "Bottom")
  openxlsx::write.xlsx(list_of_datasets, file = outexcel, asTable = TRUE, headerStyle = hs)
  
  cat('\n\n')
  cat("############ Job completed! ############\n")
  cat("OUTPUT FILES in folder:\n")
  cat("\tEXCEL: ", outexcel, "\n")
  cat("\tAbundanceLong: ",outmodeqcLong, "\n")
  cat("\tAbundanceWide: ",outmodelqc, "\n")
  cat("\tLog2fc Wide: ", outlog2fc, "\n")
  cat("\tLog2fc Impute: ", outlog2fc, "\n")
  cat("\tLog2fc Long: ", outlog2fclong, "\n")
  # cat("\tUnique per Condition: ", outUniqueProteinsCondition, "\n\n")
  
  if(enrich == "yesenrich"){
    cat("\tENRICHMENT files should also be out\n")
  }
  
}


selectTheOneInf <- function(a, b) {
  # a should be the log2fc  column
  if(a == "Inf"){
    sb <- gsub("(.*)(-)(.*)", "\\1", b)
  } else if (a == "-Inf") {
    sb <- gsub("(.*)(-)(.*)", "\\3", b)
  } else {
    sb <- 'both'
  }
  return(sb)
}

selectTheOneLog2fc <- function(a, b) {
  thrs <- 0
  # a should be the log2fc  column
  if(a > thrs){
    sb <- gsub("(.*)(-)(.*)", "\\1", b)
  } else if (a < -thrs) {
    sb <- gsub("(.*)(-)(.*)", "\\3", b)
  } else {
    sb <- 'NA'
  }
  return(sb)
}

loadModelQC = function (df_input, repro, specie) {
  
  # Remove empty entries
  if(any(df_input$PROTEIN == "")){ df_input <- df_input[-which(df_input$PROTEIN == ""),]}
  df_input$PROTEIN <- gsub("(^sp\\|)(.*)(\\|.*)", "\\2", df_input$PROTEIN )
  df_input$PROTEIN <- gsub("(.*)(\\|.*)", "\\1", df_input$PROTEIN )
  
  # TECHNICAL REPLICAS: if there are technical replicas, this means that we will find
  # two values for the same protein in the same bioreplica, therefore we need to 
  # aggregate first just in case:
  df_input <- aggregate(ABUNDANCE~PROTEIN+GROUP_ORIGINAL+SUBJECT_ORIGINAL, data = df_input, FUN = mean)
  
  dc_input <- dcast(data=df_input[,c('PROTEIN','ABUNDANCE', 'GROUP_ORIGINAL')], PROTEIN~GROUP_ORIGINAL, value.var = 'ABUNDANCE', fun.aggregate = mean, fill = 0)
  names(dc_input)[grep('PROTEIN', names(dc_input))] <- 'Protein'
  
  colnames(repro) <- paste("NumBR", colnames(repro), sep = "_")
  colnames(repro)[1]<- 'Protein'
  dc_input <- merge(dc_input, repro, by = c('Protein'))
  
  return(dc_input)
}

# Required for PCA analysis
loadModelQCstrict = function (df_input, specie) {
  
  cat("\n\t>>>Loading abundance values for proteins found in all biological replicas\n")  
  
  # Remove empty entries
  if(any(df_input$PROTEIN == "")){ df_input <- df_input[-which(df_input$PROTEIN == ""),]}
  df_input$PROTEIN <- gsub("(^sp\\|)(.*)(\\|.*)", "\\2", df_input$PROTEIN )
  df_input$PROTEIN <- gsub("(.*)(\\|.*)", "\\1", df_input$PROTEIN )
  
  # Technical replicas: aggregate on the mean the technical replicas
  b <- aggregate(ABUNDANCE~PROTEIN+GROUP_ORIGINAL+SUBJECT_ORIGINAL, data = df_input, FUN = mean)
  
  # Critical issue: BIOLOGICAL REPLICAS
  # allBiologicalReplicas <- function(x){ifelse(sum(!is.na(x)) == 3, mean(x, na.rm = T), NA)}
  
  # Before
  # datadc <- dcast(data=b, PROTEIN~GROUP_ORIGINAL, value.var = 'ABUNDANCE', fun.aggregate = mean, fill = 0)
  # before <- dim(datadc)[1]
  # l <- dim(datadc)[2]
  # datadc <- datadc[apply(datadc[c(2:l)],1,function(z) !any(z==0)),] 
  # evenafter <- dim(datadc)[1]
  # datadc <- datadc[complete.cases(datadc),]
  # after <- dim(datadc)[1]
  # cat("\t\tTotal proteins before: ", before, "\n\t\tAfter removing the 0s: ",evenafter, "\n\t\tTotal proteins (only complete cases): ", after, "\n\n")
  
  datadc <- dcast(data=b, PROTEIN~GROUP_ORIGINAL, value.var = 'ABUNDANCE', fun.aggregate = mean)  
  
  names(datadc)[grep('PROTEIN', names(datadc))] <- 'Protein'
  
  send_back <- artms_annotationUniprot(datadc, 'Protein', specie)
  return(send_back)
}

loadL2FCWide = function (df_input, repro, specie) {
  
  # # Remove the weird empty proteins
  # if(any(df_input$Protein == "")){ df_input <- df_input[-which(df_input$Protein == ""),]}
  # df_input$Protein <- gsub("(^sp\\|)(.*)(\\|.*)", "\\2", df_input$Protein )
  # df_input$Protein <- gsub("(.*)(\\|.*)", "\\1", df_input$Protein )
  
  input_melt = melt(data = df_input[,c('Protein', 'Label','log2FC','adj.pvalue'),],id.vars=c('Protein', 'Label'))
  input_dcast = dcast( Protein ~ Label+variable, data=input_melt, value.var=c('value'))
  
  colnames(repro) <- paste("NumBR", colnames(repro), sep = "_")
  colnames(repro)[1]<- 'Protein'
  input_dcast <- merge(input_dcast, repro, by = c('Protein'))
  
  # Move Gene name to the left:
  return(input_dcast)
}

loadL2FCLong = function (df_input, specie) {
  
  # Remove the weird empty proteins
  if(any(df_input$Protein == "")){ df_input <- df_input[-which(df_input$Protein == ""),]}
  df_input$Protein <- gsub("(^sp\\|)(.*)(\\|.*)", "\\2", df_input$Protein )
  df_input$Protein <- gsub("(.*)(\\|.*)", "\\1", df_input$Protein )
  
  df_input <- df_input[!is.na(df_input$log2FC),]
  df_input <- df_input[!is.infinite(df_input$log2FC),]
  df_input <- df_input[complete.cases(df_input$log2FC),]
  
  df_input$CMA <- mapply(selectTheOneLog2fc, df_input$log2FC, df_input$Label)
  
  names(df_input)[grep('Label', names(df_input))] <- 'Comparison'
  df_input$Observation <- with(df_input, paste(CMA,Comparison, sep = "--"))
  
  send_back <- df_input[,c('Protein', 'Comparison','log2FC','adj.pvalue','CMA', 'Observation')]
  
  return(send_back)
}

# POTENTIAL FUNCTION TO LOAD the ModelQC data
artms_loadModelqcBasic <- function(data){
  if( length(grep(";",data$PROTEIN))>0 ) data <- data[-grep(";",data$PROTEIN),] # NOTE!!! We lose a lot of entries this way.
  if("PROTEIN" %in% colnames(data)){
    names(data)[grep("PROTEIN", names(data))] <- 'Protein'
  }else{
    cat("ERROR: you should check the abundance file because something is seriously wrong!\n") 
    stop("Abort mission\n!")
  }
  if("ABUNDANCE" %in% colnames(data)){
    names(data)[grep("ABUNDANCE", names(data))] <- 'Abundance'
  }else{
    cat("ERROR: you should check the abundance file because something is seriously wrong!\n") 
    stop("Abort mission\n!")
  }
  if("GROUP_ORIGINAL" %in% colnames(data)){
    names(data)[grep("GROUP_ORIGINAL", names(data))] <- 'Condition'
  }else{
    cat("ERROR: you should check the abundance file because something is seriously wrong!\n") 
    stop("Abort mission\n!")
  }
  if("SUBJECT_ORIGINAL" %in% colnames(data)){
    names(data)[grep("SUBJECT_ORIGINAL", names(data))] <- 'Bioreplica'
  }else{
    cat("ERROR: you should check the abundance file because something is seriously wrong!\n") 
    stop("Abort mission\n!")
  }
  data <- subset(data, select = c(Protein, Abundance, Condition, Bioreplica))
  return(data)
}

enrichImputedLog2fc <- function(data, filenames, specie, background){
  
  # DEBUG
  # data <- filallsig_log2fc_long
  # filenames <- out.mac.allsig
  # background <- listOfGenes
  # 
  # Making sure we have unique genes in each comparison (the PTM might bring redundancy)
  pretmp <- data[c('Gene', 'Comparisons')]
  pretmp <- unique(pretmp)
  
  tmp = split(pretmp$Gene, pretmp$Comparisons, drop=T)
  
  if(specie == "human"){
    enrichgenes <- enrichProfiler(tmp, categorySource = c('GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC', 'CORUM', 'HPA','OMIM'), specie = 'hsapiens', background) # 'HP'
  }else if(specie == "mouse"){
    enrichgenes <- enrichProfiler(tmp, categorySource = c('GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC', 'CORUM'), specie = 'mmusculus', background)
  }else{
    stop("\n\n\ntOhhh no, this specie",specie," is not supported in the enrichment!!\n\n\n")
  }
  
  enrichgenes2plot <- enrichgenes[which(enrichgenes$term.size < 500),]
  # DEBUGGIN: write.table(enrichgenes2plot, 'takealookatthisshit.txt', col.names = T, row.names = F, sep = "\t", quote = F)
  for( i in unique(enrichgenes2plot$domain) ){
    cat("\t\tPlotting '", i, "' annotations...\n")
    tmp <- enrichgenes2plot[which(enrichgenes2plot$domain == i),]
    outfile <- gsub(".txt",paste0("_",i, ".txt"), filenames )
    Enrichment.plotHeatmaps(tmp, outfile)
  }
  # OUTenrich <- cleanGPROFILER(enrichgenes)
  # return(OUTenrich)
  return(enrichgenes)
}

enrichForComplexes <- function(df, backgroundNumber){
  
  listOfConditions <- unique(df$Comparisons)
  
  # fileCorum <- '~/Box Sync/db/proteinComplexes/corumKrogan20161209.txt'
  # fileCorum <- '~/Box Sync/db/proteinComplexes/coreComplexes_20170118.txt'
  fileCorum <- '~/Box Sync/db/proteinComplexes/20170801_corum_mitoT.txt'
  corumKrogan <- read.delim(file = fileCorum, header = T, sep = "\t", quote = "", stringsAsFactors = F)
  complexEnrichmentConditions <- NULL
  
  for (i in 1:length(listOfConditions)){
    condi <- listOfConditions[i]
    tmp <- unique(df$Protein[which(df$Comparisons == condi)])
    tmpEnrich <- foldComplexEnrichment(tmp, corumKrogan, backgroundNumber)
    # Check point
    checkpc <- dim(tmpEnrich)[1]
    if (checkpc > 0){
      tmpEnrich$Comparisons <- condi
      complexEnrichmentConditions <- rbind(complexEnrichmentConditions,tmpEnrich)
    }
  }
  
  return(complexEnrichmentConditions)
}

plot.corum <- function(df, outfile, theTitle){
  checkPoint <- length(unique(df$Comparisons))
  if(checkPoint >= 1){
    
    # Some of the p-values are going to be very small
    # Transform them to the smallest p-value / 10 would keep them and move it to the top
    dftemp <- df[-which(df$pvalue == 0),]
    if(dim(dftemp)[1] > 0){
      theMinimal <- min(dftemp$pvalue)/10
      # Make the value replacement
      df$p_value <- df$pvalue
      df$p_value[df$p_value == 0] <- theMinimal
    }else{
      df$p_value <- df$pvalue
    }
    
    df$p_value <- -log10(df$p_value)
    
    toplot <- dcast(data=df, ComplexName~Comparisons, value.var = "p_value", fun.aggregate = sum, fill = 0)
    rownames(toplot) <- toplot$ComplexName
    toplotmatrix <- subset(toplot, select=-c(ComplexName))
    x <- data.matrix(toplotmatrix)
    # HEATMAP
    palette.breaks <- seq(1, 3, 0.1)
    color.palette  <- colorRampPalette(c("white","steelblue"))(length(palette.breaks))
    pheatmap(x, 
             filename = outfile,
             cluster_rows = T,
             cluster_cols = F,
             cellheight = 10, 
             cellwidth=25, 
             main = theTitle,
             fontsize=6, 
             fontsize_row=8, 
             fontsize_col=12, 
             border_color='black',
             fontfamily="Helvetica",
             treeheight_row = F, 
             treeheight_col = F,
             color = color.palette
    )
    cat("\nComplex Enrichment Heatmap ready\n")
  }else{
    cat("\n!!!!!Not enough enriched comparisons to plot the heatmap\n")
  }
} #plot.corum

# IMPUTATION METHODS
# Imputing log2fc values
imputingOnTheMaximum <- function(idflog2fc){
  idflog2fc$iLog2FC <- idflog2fc$log2FC
  idflog2fc$iPvalue <- idflog2fc$adj.pvalue
  
  # this will be the loop for each condition
  theComparisons2check <- unique(idflog2fc$Label)
  
  # Getting the maximun and minimum value for each comparision
  theMaxPositive <- aggregate(data = idflog2fc[is.finite(idflog2fc$log2FC),], log2FC~Label, max)
  theMinNegative <- aggregate(data = idflog2fc[is.finite(idflog2fc$log2FC),], log2FC~Label, min)
  
  for (oneComparison in (theComparisons2check)){
    cat("\t\t Imputing", oneComparison, "\n")
    # Imputing log2fc
    idflog2fc$iLog2FC[which(idflog2fc$log2FC == "Inf" & idflog2fc$Label == oneComparison)] <- theMaxPositive$log2FC[which(theMaxPositive$Label == oneComparison)]
    idflog2fc$iLog2FC[which(idflog2fc$log2FC == "-Inf" & idflog2fc$Label == oneComparison)] <- theMinNegative$log2FC[which(theMinNegative$Label == oneComparison)]
    # Change label
    idflog2fc$imputed[which(idflog2fc$log2FC == "Inf" & idflog2fc$Label == oneComparison)] <- 'imputed'
    idflog2fc$imputed[which(idflog2fc$log2FC == "-Inf" & idflog2fc$Label == oneComparison)] <- 'imputed'
    # Imputing pvalues: for now, let's assign this minimal p-value to all the imputed values
    idflog2fc$iPvalue[which(idflog2fc$log2FC == "Inf" & idflog2fc$Label == oneComparison)] <- (0.05/numberBioReplicas)
    idflog2fc$iPvalue[which(idflog2fc$log2FC == "-Inf" & idflog2fc$Label == oneComparison)] <- (0.05/numberBioReplicas)
  }
  return(idflog2fc)
}

# NEW IMPUTATION METHOD
# When a value is completely missed in one of the conditions, 
# the log2fc = Inf / -Inf. Here, we impute those values.
# The imputation method works as follow. The assumption is that those
# proteins are likely present as well in those conditions where are missed, 
# but due to the small sampling (usually 2 or 3 biological replicas) 
# and other proteomics related issue, those proteins didn't make it through 
# the level of detection.
# Therefore, a small intensity (sampled from the bottom 5%) will be assigned
# to the protein/site in the missing condition, and the new log2fc is 
# re-calculated out of the MSstats box. Two issues are addressed in this way
# 1. If a protein has been consistently identified in one of the conditions, 
# it will stay
# 2. But if the intensity value in those conditions was too low, 
# then the log2fc will be also low

imputeMissingValues <- function(dflog2fcinfinites, dfmq) {
  
  # The comparsions
  contrast <- unique(dflog2fcinfinites$Label)
  
  # Select the IDs to impute
  ids2impute <- unique(dflog2fcinfinites$Protein)
  
  # Take the abundance values for all the proteins
  abu2imp <- artms_loadModelqcBasic(dfmq)
  # Aggregate the technical replica by choosing the maximum value
  abu2imp2 <- aggregate(Abundance~Protein+Condition+Bioreplica, data = abu2imp, FUN = mean)
  
  # Check things that will be imputed
  # dfdc.ni <- dcast(data=abu2imp2, Protein~Bioreplica, value.var = "Abundance")  
  
  # Two possible options here. 
  # 1. Select based on the bottom x%
  # # Imputing the missing values by selecting randomly from the bottom 5%
  # theMin <- min(dfmq$ABUNDANCE)
  # # Select the 5% quartile as the maximum value to sample from
  # theMax <- quantile(dfmq$ABUNDANCE, probs = .05)
  # 
  # 2. Select the bottom 20 intensities
  # Grab the bottom 30 intensities in the dataset
  dfmqOrdered <- dfmq[order(dfmq$ABUNDANCE, decreasing = F),]
  
  numberFromBottom <- 10
  abuBottom <- head(dfmqOrdered$ABUNDANCE, n = numberFromBottom)
  
  theMin <- abuBottom[1]
  theMax <- abuBottom[numberFromBottom]
  
  # Generating the numbers from which we are going to sample
  numbers2sample <- seq(from=theMin, to=theMax, by=.00001)
  
  # dcast on abundance and fill with random numbers between the minimum and q05
  suppressWarnings(dfdc.im <- dcast(data=abu2imp2, Protein~Bioreplica, value.var = "Abundance", fill = sample(numbers2sample, replace = F)  ) )
  
  # Needs to aggregate on biological replicas
  # 1. Melt on biological replicas
  dfdc.melt <- melt(dfdc.im, id.vars = c('Protein'), value.name = 'Abundance', variable.name = 'Bioreplica')
  # 2. Get the condition
  dfdc.melt$Condition <- gsub("(.*)(-)(.*)","\\1",dfdc.melt$Bioreplica)
  # 3. Dcast and Aggregate on the condition, taking the mean
  dfdc.final <- dcast(data=dfdc.melt, Protein~Condition, value.var = 'Abundance', fun.aggregate = mean)
  # 4. Filter by proteins to impute
  dfdc.final <- dfdc.final[which(dfdc.final$Protein %in% ids2impute),]
  
  for(c in contrast){
    
    cat("\t",c," --> ")
    
    x <- gsub("(.*)(-)(.*)", "\\1", c)
    y <- gsub("(.*)(-)(.*)", "\\3", c)
    cat("log2fc(",x, " - ", y,")\n")
    
    # Renaming the comparision name just for illustration purposes
    rnc <- paste0("l2fc_",c)
    
    dfdc.final[[rnc]] <- dfdc.final[[x]] - dfdc.final[[y]]
  }
  
  # Select only log2fc columns
  imputedL2FCValues <- dfdc.final[grepl("Protein|l2fc_", colnames(dfdc.final))]
  
  # Melt again
  imputedL2FCmelted <- melt(imputedL2FCValues, id.vars = c('Protein'), variable.name = 'Label', value.name = 'iLog2FC')
  # Now let's get it ready for merging with the values to be imputed at dflog2fcinfinites
  imputedL2FCmelted$Label <- gsub("l2fc_","",imputedL2FCmelted$Label)
  
  # And let's add p-values
  samplingPvalue <- seq(from=0.01, to=0.05, by=.0000001)
  # And add imputed pvalues
  imputedL2FCmelted$iPvalue <- sample(samplingPvalue, size = nrow(imputedL2FCmelted), replace = F)
  
  return (imputedL2FCmelted)
}

plotNumberProteinsImputedLog2fc <- function(data) {
  library(ggplot2)
  x <- data[c('Protein','Comparison')]
  y <- unique(x)
  z <- ggplot(y, aes(x = Comparison, fill = Comparison))
  z <- z + geom_bar(stat = "count")
  z <- z + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
  z <- z + geom_text(stat='count', aes(label=..count..), vjust=-0.5, size = 2.7)
  z <- z + ggtitle("Unique Proteins in Comparisons")  
  print(z)
}





