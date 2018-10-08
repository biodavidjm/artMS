# ------------------------------------------------------------------------------
#' @title Analysis of the Relative Quantification
#'
#' @description Analysis of the Relative Quantifications obtained by MSstats.
#' It includes:
#'
#' - Annotations
#' - Summary files in different format (xls, txt) and shapes (long, wide)
#' - Numerous summary plots
#' - Enrichment analysis using Gprofiler
#' - PCA of protein abundance
#' - PCA of quantifications
#' - Clustering analysis
#'
#' @param log2fc_file (char) MSstats results file location
#' @param modelqc_file (char) MSstats modelqc file location
#' @param specie (char) Main specie. Currently supported: human, mouse
#' @param output_dir (char) Name for the folder to output the results from the 
#' function
#' @param enrich (logical) Performed enrichment analysis using GprofileR?
#' `TRUE` (default) or `FALSE`
#' @param l2fc_thres (int) log2fc cutoff for enrichment analysis (default,
#' `l2fc_thres = 1.5`)
#' @param choosePvalue (char) specify whether `pvalue` or `adjpvalue` should 
#' use for the analysis. The default option is `adjpvalue`
#' (multiple testing correction).
#' But if the number of biological replicates for a given experiment is
#' too low (for example n = 2), then `choosePvalue = pvalue` is recommended.
#' @param isBackground (char) background of gene names for enrichment analysis.
#' `nobackground` (default) will use the total number of genes detected.
#' Alternatively provided the file path name to the background gene list.
#' @param isPtm (char) Is a ptm-site quantification? 
#' - `global` (default), 
#' - `ptmsites` (for site specific analysis), 
#' - `ptmph` (Jeff's script output evidence file)
#' @param mnbr (int) minimal number of biological replicates for imputation
#' and filtering. Default: `mnbr = 2` (Proteins must be found in one of the
#' conditions in at least 2 of the biological replicates)
#' @param isFluomics (logical) Does this data belong to the FluOMICs project?
#' `TRUE` or `FALSE` (default)
#' @param pathogen (char) Is there a pathogen in the dataset as well?
#' if it does not, then use `pathogen = nopathogen` (default).
#' Pathogens available: `tb` (Tuberculosis), `lpn` (Legionella)
#' @return (data.frame) summary of quantifications, including annotations, 
#' enrichments, etc
#' @keywords analysis, quantifications
#' @examples
#' # Testing that the files cannot be empty
#' artms_analysisQuantifications(log2fc_file = NULL,
#'                               modelqc_file = NULL,
#'                               specie = NULL,
#'                               output_dir = NULL)
#' @export
artms_analysisQuantifications <- function(log2fc_file,
                                          modelqc_file,
                                          specie,
                                          output_dir,
                                          enrich = TRUE,
                                          l2fc_thres = 1.5,
                                          choosePvalue = "adjpvalue",
                                          isBackground = "nobackground",
                                          isPtm = "global",
                                          mnbr = 2,
                                          isFluomics = FALSE,
                                          pathogen = "nopathogen") {
  cat(">> ANALYSIS OF QUANTIFICATIONS\n")
  
  if(is.null(log2fc_file) & is.null(modelqc_file) & 
     is.null(specie) & is.null(output_dir)){
    return("The evidence_file, modelqc_file, specie and output_dir arguments
must not be empty")
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Checking arguments

  # CHECK POINT: DO THE FILES EXIST?
  if(!file.exists(log2fc_file)){
    stop("THE FILE ", log2fc_file, " DOES NOT EXIST!\n")
  }
  
  if(!file.exists(modelqc_file)){
    stop("THE FILE ", modelqc_file, " DOES NOT EXIST!\n")
  }
  
  if (!grepl("logical", class(enrich))) {
    stop("\nArgument <enrich> must be logical (TRUE or FALSE)\n")
  }

  if (!grepl("logical", class(isFluomics))) {
    stop("\nArgument <isFluomics> must be logical (TRUE or FALSE)\n")
  }

  if (!grepl("numeric", class(l2fc_thres))) {
    stop("\nArgument <l2fc_thres> must be numeric\n")
  }
  
  if (!grepl("numeric", class(mnbr))) {
    stop("\nArgument <mnbr> must be numeric\n")
  }

  if(!(isPtm %in% c('global', 'ptmph', 'ptmsites'))){
    stop("The < isPtm > argument is wrong. 
         The valid options are: global or ptmsites\n")
  }
  
  if(!(choosePvalue %in% c('pvalue', 'adjpvalue'))){
    stop("The < choosePvalue > argument is wrong. 
         The valid options are: pvalue or adjpvalue\n")
  }
  
  specie <- tolower(specie)
  if(!(specie %in% c('human', 'mouse'))){
    stop("The < specie > argument is wrong. 
         The valid options are: pvalue or adjpvalue\n")
  }
  
  if (pathogen == "nopathogen") {
    cat("--- No Pathogen extra in these samples (only Influenza)\n")
  } else if (pathogen == "tb") {
    # This should not work
    cat("--- PATHOGEN IN SAMPLES: TB\n")
    pathogen.ids <- artms_data_pathogen_TB
    names(pathogen.ids) <- c('Entry')
  } else if (pathogen == "lpn") {
    cat("--- PATHOGEN IN SAMPLES: LEGIONELLA PNEUMOPHILA\n")
    pathogen.ids <- artms_data_pathogen_LPN
  } else{
    stop("\n\nThis pathogen is not supported yet\n\n")
  }
  
  output_dir <- paste0(output_dir, "_", choosePvalue)
  
  # create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # LOADING ABUNDANCE
  cat(">> LOADING modelqc FILE (ABUNDANCE)\n")
  dfmq <-
    read.delim(
      modelqc_file,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE
    )
  # Removing the empty protein names
  if (any(dfmq$PROTEIN == "")) {
    dfmq <- dfmq[-which(dfmq$PROTEIN == ""), ]
  }
  dfmq$PROTEIN <- gsub("(sp\\|)(.*)(\\|.*)", "\\2", dfmq$PROTEIN)
  dfmq$PROTEIN <- gsub("(.*)(\\|.*)", "\\1", dfmq$PROTEIN)
  # Remove outliers
  cat("--- Removing outliers (10 < abundance < 40)\n")
  dfmq <- dfmq[which(dfmq$ABUNDANCE > 10 & dfmq$ABUNDANCE < 40), ]
  
  # First, let's take the conditions, which will be used later in several places
  conditions <- unique(dfmq$GROUP_ORIGINAL)
  numberConditions <- length(conditions)
  
  # KEY STEP: GETTING THE BACKGROUND GENE LIST
  if (isBackground == "nobackground") {
    # If not list of background genes is provided,
    # then extract them from the modelqc file
    if (isPtm == "global") {
      suppressMessages(dfmq2Genes <-
                         artms_annotationUniprot(dfmq, 'PROTEIN', specie))
      numberTotalGenes <- length(unique(dfmq2Genes$Gene))
      cat("--- TOTAL NUMBER OF GENES/PROTEINS: ",
          numberTotalGenes,
          "\n")
      listOfGenes <- unique(dfmq2Genes$Gene)
    } else if (grepl("ptm", isPtm)) {
      dfmq2Genes <-
        dfmq[c('PROTEIN', 'GROUP_ORIGINAL')] 
      names(dfmq2Genes)[grep('PROTEIN', names(dfmq2Genes))] <- 'Protein'
      # Removing party sites
      dfmq2Genes <-
        dfmq2Genes[grep(",", dfmq2Genes$Protein, invert = TRUE), ]
      # And now be very careful with the Fluomics labeling, since they have an 
      # extra _ that it is not follow by the site
      cat(
        "--- Warning! if you have UNIPROT_PTM id with more than one 
        underscore '_' is going to be a problem\n"
      )
      dfmq2Genes$Protein <- ifelse(
        grepl("_H1N1|_H3N2|_H5N1", dfmq2Genes$Protein),
        gsub(
          "^(\\S+?_H[1,3,5]N[1,2])_.*",
          "\\1",
          dfmq2Genes$Protein,
          perl = TRUE
        ) ,
        gsub("^(\\S+?)_.*", "\\1", dfmq2Genes$Protein, perl = TRUE)
      )
      dfmq2Genes <- unique(dfmq2Genes)
      suppressMessages(dfmq2Genes <-
                         artms_annotationUniprot(dfmq2Genes, 'Protein', specie))
      numberTotalGenes <- length(unique(dfmq2Genes$Gene))
      cat("--- TOTAL NUMBER OF GENES/PROTEINS: ",
          numberTotalGenes,
          "\n\n")
      if (numberTotalGenes == 0) {
        stop("\nSOMETHING WRONG WITH THE IDs. Check the source code\n")
      }
      listOfGenes <- unique(dfmq2Genes$Gene)
    }
  } else{
    # No matter what list is provided, it must come with a "Gene" column
    backgroundList <-
      read.delim(
        isBackground,
        header = TRUE,
        sep = "\t",
        quote = "",
        stringsAsFactors = FALSE
      )
    listOfGenes <- unique(backgroundList$Gene)
  }
  
  backgroundNumber <- length(listOfGenes)
  
  #-----------------------------------------------------------------------------
  # LOG2FC
  cat(">> LOADING QUANTIFICATIONS (-results.txt from MSstats)\n")
  dflog2fcraw <-
    read.delim(
      log2fc_file,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE
    )
  if (any(dflog2fcraw$Protein == "")) {
    dflog2fcraw <- dflog2fcraw[-which(dflog2fcraw$Protein == ""), ]
  }
  dflog2fcraw$Protein <-
    gsub("(sp\\|)(.*)(\\|.*)", "\\2", dflog2fcraw$Protein)
  dflog2fcraw$Protein <-
    gsub("(.*)(\\|.*)", "\\1", dflog2fcraw$Protein)
  
  # # Filtering conditions:
  # if(isFluomics == TRUE){
  #   cat("WARNING! selecting MOCK for humans and mice in LOG2FC raw data\n")
  #
  #   # WHEN A REFERENCE MOCK IS WISHED
  #   # if(specie == "human"){
  #   #   dflog2fcraw <- dflog2fcraw[(grepl("H[[:digit:]]N[[:digit:]]", 
  #   dflog2fcraw$Label) & grepl("MOCK_03H", dflog2fcraw$Label) ),]
  #   # }else if (specie == "mouse"){
  #   #   dflog2fcraw <- dflog2fcraw[(grepl("H[[:digit:]]N[[:digit:]]", 
  #   dflog2fcraw$Label) | grepl("MOCK_D04", dflog2fcraw$Label) ),]
  #   # }
  #
  #   #Choosing the match mock
  #   flu_contrast <- c("H1N1_03H-MOCK_03H", "H1N1_06H-MOCK_06H", 
  #   "H1N1_12H-MOCK_12H", "H1N1_18H-MOCK_18H", "H3N2_03H-MOCK_03H", 
  #   "H3N2_06H-MOCK_06H", "H3N2_12H-MOCK_12H", "H3N2_18H-MOCK_18H", 
  #   "H5N1_03H-MOCK_03H", "H5N1_06H-MOCK_06H", "H5N1_12H-MOCK_12H", 
  #   "H5N1_18H-MOCK_18H")
  #   dflog2fcraw <- dflog2fcraw[which(dflog2fcraw$Label %in% flu_contrast),]
  #
  #   cat("LOG2FC Data Filtered by specific FLU comparisons\n")
  # }
  
  # Let's get rid of outliers: log2fc larger than X (but we need to keep 
  # the "inf" values for imputation)
  dflog2fcfinites <- dflog2fcraw[is.finite(dflog2fcraw$log2FC), ]
  cutofflog2fc <- 14
  cat(
    "--- Removing outliers (",
    paste0("-", cutofflog2fc, " < log2fc < +", cutofflog2fc),
    ")\n"
  )
  filtermorethan10 <-
    length(dflog2fcfinites$log2FC[abs(dflog2fcfinites$log2FC) > cutofflog2fc])
  if (filtermorethan10 > 0) {
    cat(
      "------ (-) Removing [",
      filtermorethan10,
      "] protein ids with a abs(log2fc) >",
      cutofflog2fc,
      "\n"
    )
    dflog2fcfinites <-
      dflog2fcfinites[-which(abs(dflog2fcfinites$log2FC) > cutofflog2fc), ]
  } else{
    cat("------ (+) No abs(log2fc) >", cutofflog2fc, "so moving on!\n")
  }
  
  # IMPUTING MISSING VALUES
  # When a value is completely missed in one of the conditions,
  # the log2fc = Inf / -Inf. Here, we impute those values.
  # The imputation method works as follow. The assumption is that those
  # proteins are likely present as well in those conditions where are missed, 
  # but due to the
  # small sampling (usually 2 or 3 biological replicas) and other proteomics
  # related issue, those proteins didn't make it through the level of detection.
  # Therefore, a small intensity (sampled from the bottom 5%) will be assigned
  # to the protein/site in the missing condition, and the new log2fc is 
  # re-calculated out of the MSstats box. Two issues are addressed in this way
  # 1. If a protein has been consistently identified in one of the conditions, 
  # it will stay
  # 2. But if the intensity value in those conditions was too low, then the 
  # log2fc will be also low
  
  cat(">> IMPUTING MISSING VALUES\n")
  # Select infinite values (i.e., log2fc missed for that)
  dflog2fcinfinites <- dflog2fcraw[is.infinite(dflog2fcraw$log2FC), ]
  numberInfinites <- dim(dflog2fcinfinites)[1]
  
  # Control
  if (numberInfinites < 1) {
    cat("\nWARNING: O infinite values. This is not normal\n")
  } else {
    cat("--- Number of +/- INF values", dim(dflog2fcinfinites)[1], "\n")
    
    imputedL2FCmelted <-
      .artms_imputeMissingValues(dflog2fcinfinites, dfmq)
    
    # Merge with the original log2fc values to impute...
    theImputedL2FC <-
      merge(
        dflog2fcinfinites,
        imputedL2FCmelted,
        by = c("Protein", "Label"),
        all.x = TRUE
      )
    theImputedL2FC$imputed <- "yes"
  }
  
  # Getting the data ready for merging
  dflog2fcfinites$imputed <- "no"
  dflog2fcfinites$iLog2FC <- dflog2fcfinites$log2FC
  
  # Choose the pvalue or adjusted pvalue as the iPvalue
  if (choosePvalue == "pvalue") {
    dflog2fcfinites$iPvalue <- dflog2fcfinites$pvalue
  } else if (choosePvalue == "adjpvalue") {
    dflog2fcfinites$iPvalue <- dflog2fcfinites$adj.pvalue
  } else{
    stop("\n\n\t------> wait a minute: did you choose pvalue or adjpvalue")
  }
  
  # Merging: NA values are thrown away at this point
  if (numberInfinites < 1) {
    dflog2fc <- dflog2fcfinites
  } else {
    dflog2fc <- rbind(dflog2fcfinites, theImputedL2FC)
  }
  
  cat("--- Plotting distributions of log2fc and pvalues\n")
  
  plotDFdistColor <-
    ggplot(dflog2fc, aes(x = log2FC, fill = Label)) +
    geom_histogram(bins = 100,
                   alpha = .4,
                   col = "black") +
    labs(title = "Distribution log2FC", x = "log2FC")
  
  plotDFdistAll <- ggplot(dflog2fc, aes(x = log2FC)) +
    geom_histogram(bins = 100,
                   alpha = .4,
                   col = "black") +
    labs(title = "Distribution log2FC", x = "log2FC")
  
  plotDFdistiLog <- ggplot(dflog2fc, aes(x = iLog2FC)) +
    geom_histogram(bins = 100,
                   alpha = .4,
                   col = "black") +
    labs(title = "Distribution ilog2FC (imputed + nonimputed", x = "iLog2FC")
  
  plotPvalues <-
    ggplot(dflog2fc[is.finite(dflog2fc$pvalue), ], aes(x = pvalue)) +
    geom_histogram(bins = 50,
                   alpha = .4,
                   col = "black") +
    labs(title = "Distribution p-values", x = "p-values")
  
  plotAdjustedPvalues <-
    ggplot(dflog2fc[-which(dflog2fc$adj.pvalue == 0), ], aes(x = adj.pvalue)) +
    geom_histogram(bins = 150,
                   alpha = .4,
                   col = "black") +
    labs(title = "Distribution adj.pvalues", x = "adj.values")
  
  plotAdjustedIpvalues <- ggplot(dflog2fc, aes(x = iPvalue)) +
    geom_histogram(bins = 150,
                   alpha = .4,
                   col = "black") +
    labs(title = "Distribution imputed p-values", x = "iPvalues")
  
  
  # DISTRIBUTION PRINT OUTS
  distributionsOut <- gsub(".txt", ".distributions.pdf", log2fc_file)
  distributionsOut <- paste0(output_dir, "/", distributionsOut)
  pdf(distributionsOut)
  plotDFdistColor
  plotDFdistAll
  plotDFdistiLog
  plotPvalues
  plotAdjustedPvalues
  plotAdjustedIpvalues
  
  if (numberInfinites > 0) {
    hist(
      imputedL2FCmelted$iLog2FC,
      breaks = 100,
      main = paste0("Imputed Log2FC (all)\n  n = ", dim(imputedL2FCmelted)[1]),
      xlab = "log2fc"
    )
    hist(
      theImputedL2FC$iLog2FC,
      breaks = 100,
      main = paste0("Imputed Log2FC merged\n n = ", dim(theImputedL2FC)[1]),
      xlab = "log2fc"
    )
  }
  
  hist(
    dflog2fcfinites$pvalue,
    breaks = 100,
    main = paste0("p-value distribution\n n = ", dim(dflog2fcfinites)[1]),
    xlab = "adj.pvalues"
  )
  hist(
    dflog2fcfinites$adj.pvalue,
    breaks = 100,
    main = paste0("Adjusted p-values distribution\n n = ", 
                  dim(dflog2fcfinites)[1]),
    xlab = "adj.pvalues"
  )
  hist(
    dflog2fcfinites$iLog2FC,
    breaks = 1000,
    main = paste0(
      "Non-imputed Log2FC distribution\n n = ",
      dim(dflog2fcfinites)[1]
    ),
    xlab = "log2FC"
  )
  hist(
    dflog2fc$iPvalue,
    breaks = 100,
    main = paste0(
      "(Imputed+NonImputed) adjusted pvalue distribution\n n = ",
      dim(dflog2fc)[1]
    ),
    xlab = "adj.pvalues"
  )
  hist(
    dflog2fc$iLog2FC,
    breaks = 1000,
    main = paste0(
      "(Imputed+NonImputed) log2fc distribution\n n = ",
      dim(dflog2fc)[1]
    ),
    xlab = "log2FC"
  )
  garbage <- dev.off()
  
  
  # Relationship between conditions
  # Get the number of biological replicas based on the first condition
  theConditions <- unique(dfmq$GROUP_ORIGINAL)
  theFirstCond <- theConditions[2]
  condFirst <- dfmq[which(dfmq$GROUP_ORIGINAL == theFirstCond), ]
  theBiologicalReplicas <- unique(condFirst$SUBJECT_ORIGINAL)
  numberBioReplicas <- length(theBiologicalReplicas)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ABUNDANCE PLOTS
  
  cat(">> PLOTS: ABUNDANCE PLOTS\n")
  abundancesName <-
    gsub(".txt", ".relativeABUNDANCE.pdf", log2fc_file)
  abundancesName <- paste0("plot.", abundancesName)
  abundancesName <- paste0(output_dir, "/", abundancesName)
  pdf(abundancesName)
  .artms_plotAbundanceBoxplots(data = dfmq)
  .artms_plotNumberProteinsAbundance(data = dfmq)
  garbage <- dev.off()
  
  # Reproducibility plots based on normalized abundance
  cat(">> PLOTS: REPRODUCIBILITY PLOTS\n")
  reproName <-
    gsub(".txt", ".reproducibilityAbundance.pdf", log2fc_file)
  reproName <- paste0("plot.", reproName)
  reproName <- paste0(output_dir, "/", reproName)
  pdf(reproName)
  .artms_plotReproducibilityAbundance(dfmq)
  garbage <- dev.off()
  
  # Conditions
  cat(">> PLOT: CORRELATION BETWEEN ALL CONDITIONS\n")
  relaCond <-
    gsub(".txt", ".correlationConditions.pdf", log2fc_file)
  relaCond <- paste0("plot.", relaCond)
  relaCond <- paste0(output_dir, "/", relaCond)
  pdf(relaCond)
  .artms_plotCorrelationConditions(dfmq, numberBioReplicas)
  garbage <- dev.off()
  
  # Relationship between log2fc comparisons
  cat(">> PLOT: CORRELATION BETWEEN QUANTIFICATIONS (based on log2fc values\n")
  if (length(unique(dflog2fc$Label)) > 1) {
    relaChanges <-
      gsub(".txt", ".correlationQuantifications.pdf", log2fc_file)
    relaChanges <- paste0("plot.", relaChanges)
    relaChanges <- paste0(output_dir, "/", relaChanges)
    pdf(relaChanges)
    .artms_plotRatioLog2fc(dflog2fc)
    garbage <- dev.off()
  } else{
    cat("--- Only one Comparison is available (correlation is not possible)\n")
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ABUNDANCE DATA, CREATE FILTERS
  abundance <- .artms_loadModelqcBasic(dfmq)
  names(abundance)[grep('Protein', names(abundance))] <- 'Prey'
  names(abundance)[grep('Condition', names(abundance))] <- 'Bait'
  # TECHNICAL REPLICAS: if there are technical replicas means that we will find
  # two values for the same protein in the same bioreplica, therefore we need to
  # aggregate first just in case:
  abundance <-
    aggregate(Abundance ~ Prey + Bait + BioReplicate,
              data = abundance,
              FUN = mean)
  
  # Let's aggregate to get the sum of the abundance, we will use it later.
  abundance_dcsum <-
    data.table::dcast(
      abundance,
      Prey ~ Bait,
      value.var = 'Abundance',
      fun.aggregate = sum,
      fill = 0
    )
  abundance_dcmean <-
    data.table::dcast(
      abundance,
      Prey ~ Bait,
      value.var = 'Abundance',
      fun.aggregate = mean,
      fill = 0
    )
  
  #########################################################
  # HEATMAPs: Use the mean for the heatmap
  cat(">> HEATMAPS OF PROTEIN ABUNDANCE\n")
  dchm_input <- abundance_dcmean
  rownames(dchm_input) <- dchm_input$Prey
  dfhm <- subset(dchm_input, select = -c(Prey))
  aqui <- data.matrix(dfhm)
  
  outHeatMapOverall <-
    gsub(".txt",
         ".clustering.abundance.all-overview.pdf",
         log2fc_file)
  outHeatMapOverall <- paste0(output_dir, "/", outHeatMapOverall)
  pheatmap(
    aqui,
    filename = outHeatMapOverall,
    cellwidth = 20,
    main = "Clustered Relative Abundance",
    cluster_cols = FALSE,
    fontfamily = "Helvetica",
    labels_row = "",
    fontsize = 6,
    fontsize_row = 8,
    fontsize_col = 8,
    border_color = NA,
    fontfamily = "Helvetica"
  )
  outHeatMapZoom <-
    gsub(".txt", ".clustering.abundance.all-zoom.pdf", log2fc_file)
  outHeatMapZoom <- paste0(output_dir, "/", outHeatMapZoom)
  pheatmap(
    aqui,
    filename = outHeatMapZoom,
    cellheight = 10,
    cellwidth = 20,
    main = "Clustered Relative Abundance",
    cluster_cols = FALSE,
    fontsize = 6,
    fontsize_row = 8,
    fontsize_col = 8,
    border_color = NA,
    fontfamily = "Helvetica"
  )
  
  # Melt again the sum and mean
  abundancelongsum <-
    reshape2::melt(
      abundance_dcsum,
      id.vars = c('Prey'),
      value.name = 'Abundance',
      variable.name = 'Bait'
    )
  abundancelongmean <-
    reshape2::melt(
      abundance_dcmean,
      id.vars = c('Prey'),
      value.name = 'Abundance',
      variable.name = 'Bait'
    )
  
  # We dont need the 0 values
  abundancelongsum <-
    abundancelongsum[!(abundancelongsum$Abundance == 0), ]
  abundancelongmean <-
    abundancelongmean[!(abundancelongmean$Abundance == 0), ]
  # Rename and merge:
  names(abundancelongsum)[grep('Abundance', names(abundancelongsum))] <-
    'AbSum'
  names(abundancelongmean)[grep('Abundance', names(abundancelongmean))] <-
    'AbMean'
  
  abundancelongsummean <-
    merge(abundancelongsum, abundancelongmean, by = c('Prey', 'Bait'))
  
  # REPRODUCIBILITY AND SPECIFICY PARATEMERS
  
  # Get the number of bioreplicates based on abundance data
  nbr_wide <-
    data.table::dcast(
      abundance,
      Prey ~ Bait,
      value.var = 'Abundance',
      fun.aggregate = length,
      fill = 0
    )
  nbr_long <-
    reshape2::melt(
      nbr_wide,
      id.vars = c('Prey'),
      value.name = 'Abundance',
      variable.name = 'Bait'
    )
  nbr_long <- nbr_long[!(nbr_long$Abundance == 0), ]
  names(nbr_long)[grep('Abundance', names(nbr_long))] <- 'BioRep'
  
  # Get
  OUTreprod <-
    data.table::dcast(data = nbr_long, Prey ~ Bait, value.var = 'BioRep')
  here <- dim(OUTreprod)[2]
  OUTreprod[is.na(OUTreprod)] <- 0
  # Make a copy to use later
  bioReplicaInfo <- OUTreprod
  # And geht the total number of biological replicates
  OUTreprod$ReproBioreplicaCount <- rowSums(OUTreprod[, 2:here])
  
  # Get whether a protein is found in all conditions
  reprospec2merge <-
    subset(OUTreprod, select = c(Prey, ReproBioreplicaCount))
  
  OUTreproCondition <-
    data.table::dcast(data = nbr_long, Prey ~ Bait, value.var = 'BioRep')
  thedim <- dim(OUTreproCondition)[2]
  OUTreproCondition[is.na(OUTreproCondition)] <- 0
  thepreys <- subset(OUTreproCondition, select = c(Prey))
  thevalues <- subset(OUTreproCondition, select = -c(Prey))
  thevalues[thevalues > 0] <- 1
  thevalues$ReproConditionCount <- rowSums(thevalues)
  FinalReproCondition <- cbind(thepreys, thevalues)
  
  reprocondition2merge <-
    subset(FinalReproCondition, select = c(Prey, ReproConditionCount))
  
  # This version will be printed out below
  OUTreprodFinal <-
    merge(OUTreprod, reprocondition2merge, by = 'Prey')
  
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # PCA ANALYSIS
  # It requires a simplified version for modelqc
  cat(">> PRINCIPAL COMPONENT ANALYSIS BASED ON ABUNDANCE\n")
  modelqcabundance <-
    .artms_loadModelQCstrict(df_input = dfmq,
                             specie = specie,
                             ptmis = isPtm)
  out.pca <- gsub(".txt", "-pca", log2fc_file)
  out.pca <- paste0(output_dir, "/", out.pca)
  suppressWarnings(.artms_getPCAplots(modelqcabundance, out.pca, conditions))
  cat("---+ PCA done!\n")
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ANNOTATIONS
  modelqc_file_splc <- .artms_mergeAbNbr(dfmq, nbr_wide, specie)
  
  cat(">>> ANNOTATIONS\n")
  # Now get ready for annotation
  cat("--- Abundance data\n")
  if (grepl("ptm", isPtm)) {
    names(modelqc_file_splc)[grep('^Protein$', names(modelqc_file_splc))] <-
      'Uniprot_PTM'
    # Take the Protein ID, but being very careful about the fluomics labeling
    modelqc_file_splc$Protein <-
      ifelse(
        grepl("_H1N1|_H3N2|_H5N1", modelqc_file_splc$Uniprot_PTM),
        gsub(
          "^(\\S+?_H[1,3,5]N[1,2])_.*",
          "\\1",
          modelqc_file_splc$Uniprot_PTM,
          perl = TRUE
        ) ,
        gsub(
          "^(\\S+?)_.*",
          "\\1",
          modelqc_file_splc$Uniprot_PTM,
          perl = TRUE
        )
      )
    suppressMessages(
      modelqc_file_splc <-
        artms_annotationUniprot(modelqc_file_splc, 'Protein', specie)
    )
  } else{
    suppressMessages(
      modelqc_file_splc <-
        artms_annotationUniprot(modelqc_file_splc, 'Protein', specie)
    )
  }
  
  cat("--- Relative Quantifications (Log2fc)\n")
  # Prepare output of changes
  log2fc_file_splc <-
    .artms_mergeChangesNbr(dflog2fc, nbr_wide, specie)
  # Now get ready for annotation
  if (grepl("ptm", isPtm)) {
    names(log2fc_file_splc)[grep('^Protein$', names(log2fc_file_splc))] <-
      'Uniprot_PTM'
    # Take the Protein ID, but being very careful about the fluomics labeling
    log2fc_file_splc$Protein <-
      ifelse(
        grepl("_H1N1|_H3N2|_H5N1", log2fc_file_splc$Uniprot_PTM),
        gsub(
          "^(\\S+?_H[1,3,5]N[1,2])_.*",
          "\\1",
          log2fc_file_splc$Uniprot_PTM,
          perl = TRUE
        ) ,
        gsub("^(\\S+?)_.*", "\\1", log2fc_file_splc$Uniprot_PTM, perl = TRUE)
      )
    suppressMessages(
      log2fc_file_splc <-
        artms_annotationUniprot(log2fc_file_splc, 'Protein', specie)
    )
  } else{
    suppressMessages(
      log2fc_file_splc <-
        artms_annotationUniprot(log2fc_file_splc, 'Protein', specie)
    )
  }
  
  cat(">> FILTERING CHANGES BEFORE PRINTING OUT\n")
  
  imputedDF <-
    dflog2fc[c(
      'Protein',
      'Label',
      'log2FC',
      'pvalue',
      'adj.pvalue',
      'imputed',
      'iLog2FC',
      'iPvalue'
    )]
  cat("--- Merging Changes with bioReplica Info\n")
  imputedDF <-
    merge(
      imputedDF,
      bioReplicaInfo,
      by.x = 'Protein',
      by.y = 'Prey',
      all.x = TRUE
    )
  
  cat("--- Removing NA\n")
  imputedDF <- imputedDF[!is.na(imputedDF$log2FC), ]
  
  cat("--- Add labeling of condition more abundant in the quantification\n")
  imputedDF$CMA <-
    mapply(.artms_selectTheOneLog2fc,
           imputedDF$iLog2FC,
           imputedDF$Label)
  
  cat(
    "--- Removing proteins not found in a minimal number (",
    mnbr,
    ") of biological replicates\n"
  )
  imputedDF <- .artms_RemoveProtBelowThres(imputedDF, mnbr)
  
  cat("--- Filtering is done!\n")
  
  cat(">> GENERATING QC PLOTS ABOUT CHANGES (log2fc)\n")
  cat("--- Distribution of log2fc and pvalues\n")
  distributionsFilteredOut <-
    gsub(".txt", ".distributionsFil.pdf", log2fc_file)
  distributionsFilteredOut <-
    paste0(output_dir, "/", distributionsFilteredOut)
  pdf(distributionsFilteredOut)
  hist(
    imputedDF$iLog2FC,
    breaks = 1000,
    main = paste0("Filtered Log2FC (>2BR)\n n = ", dim(imputedDF)[1]),
    xlab = "log2fc"
  )
  hist(
    imputedDF$iPvalue,
    breaks = 1000,
    main = paste0("Filtered p-values (>2BR)\n n = ", dim(imputedDF)[1]),
    xlab = "p-value"
  )
  garbage <- dev.off()
  
  cat("--- Proportion imputed values\n")
  # Stats about imputed values
  yesimputed <- dim(imputedDF[which(imputedDF$imputed == 'yes'), ])[1]
  nonimputed <- dim(imputedDF[which(imputedDF$imputed == 'no'), ])[1]
  
  dat <-
    data.frame(
      count = c(yesimputed, nonimputed),
      category = c("Imputed", "Non-Imputed")
    )
  # Add addition columns, needed for drawing with geom_rect.
  dat$fraction = dat$count / sum(dat$count)
  dat <- dat[order(dat$fraction),]
  dat$ymax <- cumsum(dat$fraction)
  dat$ymin <- c(0, head(dat$ymax, n = -1))
  
  p1 <-
    ggplot(dat, aes(
      fill = category,
      ymax = ymax,
      ymin = ymin,
      xmax = 4,
      xmin = 3
    )) +
    geom_rect() +
    coord_polar(theta = "y") +
    xlim(c(0, 4)) +
    labs(title = "Proportion of Imputed Intensity values")
  p2 <- ggplot(dat, aes(x = category, y = count, fill = category)) +
    geom_bar(stat = "identity") +
    labs(title = "Proportion of Imputed Intensity values")
  
  outImputation <- gsub(".txt", ".imputation.pdf", log2fc_file)
  outImputation <- paste0(output_dir, "/", outImputation)
  pdf(outImputation)
  print(p1)
  print(p2)
  garbage <- dev.off()
  
  cat(">> HEATMAPS OF CHANGES (log2fc)\n")
  l2fcol <-
    data.table::dcast(data = imputedDF, Protein ~ Label, value.var = 'iLog2FC')
  rownames(l2fcol) <- l2fcol$Protein
  l2fcol <- within(l2fcol, rm(Protein))
  l2fcol[is.na(l2fcol)] <- 0
  
  if (numberConditions > 1) {
    l2fcolmatrix <- data.matrix(l2fcol)
    cat("--- All changes\n")
    outHeatMapOverallL2fc <-
      gsub(".txt",
           ".clustering.log2fc.all-overview.pdf",
           log2fc_file)
    outHeatMapOverallL2fc <-
      paste0(output_dir, "/", outHeatMapOverallL2fc)
    pheatmap(
      l2fcolmatrix,
      filename = outHeatMapOverallL2fc,
      cellwidth = 20,
      main = "Clustering Log2FC",
      cluster_cols = FALSE,
      clustering_method = "average",
      fontfamily = "Helvetica",
      show_colnames = FALSE,
      fontsize = 6,
      fontsize_row = 3,
      fontsize_col = 10,
      border_color = NA
    )
    outHeatMapZoomL2fc <-
      gsub(".txt", ".clustering.log2fc.all-zoom.pdf", log2fc_file)
    outHeatMapZoomL2fc <- paste0(output_dir, "/", outHeatMapZoomL2fc)
    pheatmap(
      l2fcolmatrix,
      filename = outHeatMapZoomL2fc,
      cellheight = 10,
      cellwidth = 20,
      main = "Clustering Log2FC",
      cluster_cols = FALSE,
      fontsize = 6,
      fontsize_row = 8,
      fontsize_col = 8,
      border_color = NA,
      fontfamily = "Helvetica"
    )
    
    # Only significant pvalues
    cat("--- Only significant changes\n")
    imputedDFsig <- imputedDF[which(imputedDF$iPvalue < 0.05), ]
    l2fcolSignificants <-
      data.table::dcast(data = imputedDFsig, 
                        Protein ~ Label, 
                        value.var = 'iLog2FC')
    rownames(l2fcolSignificants) <- l2fcolSignificants$Protein
    l2fcolSignificants <- within(l2fcolSignificants, rm(Protein))
    l2fcolSignificants[is.na(l2fcolSignificants)] <- 0
    
    l2fcolSignificantsmatrix <- data.matrix(l2fcolSignificants)
    outHeatMapOverallL2fc <-
      gsub(".txt",
           ".clustering.log2fcSign.all-overview.pdf",
           log2fc_file)
    outHeatMapOverallL2fc <-
      paste0(output_dir, "/", outHeatMapOverallL2fc)
    pheatmap(
      l2fcolSignificantsmatrix,
      filename = outHeatMapOverallL2fc,
      cellwidth = 20,
      main = "Clustering Log2FC (p-value < 0.05)",
      cluster_cols = FALSE,
      fontfamily = "Helvetica",
      labels_row = "",
      fontsize = 6,
      fontsize_row = 8,
      fontsize_col = 8,
      border_color = NA,
      fontfamily = "Helvetica"
    )
    outHeatMapZoomL2fc <-
      gsub(".txt",
           ".clustering.log2fcSign.all-zoom.pdf",
           log2fc_file)
    outHeatMapZoomL2fc <- paste0(output_dir, "/", outHeatMapZoomL2fc)
    pheatmap(
      l2fcolSignificantsmatrix,
      filename = outHeatMapZoomL2fc,
      cellheight = 10,
      cellwidth = 20,
      main = "Clustering Log2FC (p-value < 0.05)",
      cluster_cols = FALSE,
      fontsize = 6,
      fontsize_row = 8,
      fontsize_col = 8,
      border_color = NA,
      fontfamily = "Helvetica"
    )
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ENRICHMENT OF MOST ABUNDANT PROTEINS (from IMPUTED LOG2FC values)
  # Let's select first significance based on pvalue, by using the iPvalue
  # we are already including the imputed pvalues...
  l2fcol4enrichment <-
    data.table::dcast(data = imputedDF[which(imputedDF$iPvalue < 0.05), ], 
                      Protein ~ Label, value.var = 'iLog2FC')
  
  if (dim(l2fcol4enrichment)[1] > 0) {
    # Let's melt now for enrichment analysis
    l2fcol4enrichment <-
      reshape2::melt(
        data = l2fcol4enrichment,
        id.vars = c('Protein'),
        variable.name = "Comparisons"
      )
    l2fcol4enrichment <-
      l2fcol4enrichment[complete.cases(l2fcol4enrichment), ]
    # Now get ready for annotation
    if (grepl("ptm", isPtm)) {
      names(l2fcol4enrichment)[grep('^Protein$', names(l2fcol4enrichment))] <-
        'Uniprot_PTM'
      # Take the Protein ID, but being very careful about the fluomics labeling
      l2fcol4enrichment$Protein <-
        ifelse(
          grepl("_H1N1|_H3N2|_H5N1", l2fcol4enrichment$Uniprot_PTM),
          gsub(
            "^(\\S+?_H[1,3,5]N[1,2])_.*",
            "\\1",
            l2fcol4enrichment$Uniprot_PTM,
            perl = TRUE
          ) ,
          gsub(
            "^(\\S+?)_.*",
            "\\1",
            l2fcol4enrichment$Uniprot_PTM,
            perl = TRUE
          )
        )
      suppressMessages(
        l2fcol4enrichment <-
          artms_annotationUniprot(l2fcol4enrichment, 'Protein', specie)
      )
    } else{
      suppressMessages(
        l2fcol4enrichment <-
          artms_annotationUniprot(l2fcol4enrichment, 'Protein', specie)
      )
    }
  }
  
  if (enrich == TRUE & dim(l2fcol4enrichment)[1] > 0) {
    cat(">> ENRICHMENT ANALYSIS OF SELECTED CHANGES USING GPROFILER\n")
    
    if (grepl("ptm", isPtm)) {
      # l2fcol4enrichment <- 
      # within(l2fcol4enrichment, rm(Gene,Uniprot_PTM,Protein.names))
      # Remove parties for enrichment
      l2fcol4enrichment <-
        l2fcol4enrichment[grep(",", l2fcol4enrichment$Protein, invert = TRUE),]
      # Select the Uniprot ID, but keep in mind that some of them might
      # have many _ph54_ph446 before
      # l2fcol4enrichment$Protein <- 
      # gsub("^(\\S+?)_.*", "\\1", l2fcol4enrichment$Protein, perl = TRUE)
      l2fcol4enrichment <-
        unique(l2fcol4enrichment[c("Protein", "Gene", "Comparisons", "value")])
    }
    
    # ALL SIGNIFICANT CHANGES log2fc
    # GPROFILER
    cat("1) Enrichment of ALL significant Changes\n")
    
    filallsig_log2fc_long <-
      l2fcol4enrichment[which(abs(l2fcol4enrichment$value) >= l2fc_thres),]
    
    if (dim(filallsig_log2fc_long)[1] > 0) {
      out.mac.allsig <-
        gsub(".txt", "-enrich-MAC-allsignificants.txt", log2fc_file)
      out.mac.allsig <- paste0(output_dir, "/", out.mac.allsig)
      mac.allsig <- NULL
      
      tryCatch(
        mac.allsig <- artms_enrichLog2fc(
          dataset = filallsig_log2fc_long,
          output_name = out.mac.allsig,
          specie = specie,
          heatmaps = TRUE,
          background = listOfGenes
        ), error = function(e){
          cat("\n\n------ (!! Error): Enrichment is not possible!\n")
          cat("                    gProfiler server is likely down\n")
          cat("                    Just wait until is up again\n\n")
          }
      )
      
      if(!is.null(mac.allsig)){
        if (dim(mac.allsig)[1] > 0) {
          write.table(
            mac.allsig,
            out.mac.allsig,
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE
          )
        }
      }

      
      cat("---+ Corum Protein Complex Enrichment Analysis\n")
      
      # CORUM
      allsigComplexEnriched <-
        .artms_enrichForComplexes(filallsig_log2fc_long, backgroundNumber)
      
      if (dim(allsigComplexEnriched)[1] > 0) {
        out.mac.allsig.corum <-
          gsub(".txt",
               "-enrich-MAC-allsignificants-corum.txt",
               log2fc_file)
        out.mac.allsig.corum <-
          paste0(output_dir, "/", out.mac.allsig.corum)
        write.table(
          allsigComplexEnriched,
          out.mac.allsig.corum,
          quote = FALSE,
          sep = "\t",
          row.names = FALSE,
          col.names = TRUE
        )
      }
      # And the heatmap
      if (dim(allsigComplexEnriched)[1] > 2) {
        out.mac.allsig.corum.pdf <-
          gsub(".txt",
               "-enrich-MAC-allsignificants-corum.pdf",
               log2fc_file)
        out.mac.allsig.corum.pdf <-
          paste0(output_dir, "/", out.mac.allsig.corum.pdf)
        .artms_plotCorumEnrichment(
          allsigComplexEnriched,
          out.mac.allsig.corum.pdf,
          "MAC ALL SIGNIFICANT Protein Complex Enrichment"
        )
      } else{
        cat("--- (-) Not enough negative corum complexes to plot\n")
      }
    } else{
      cat("\n----(-) NOTHING is significant! Check what's going on\n\n")
      mac.allsig <- NULL
      allsigComplexEnriched <- NULL
    }
    
    # POSITIVE log2fc
    # GPROFILER
    cat("2) Enrichment of selected POSITIVE significant changes\n")
    
    filpos_log2fc_long <-
      l2fcol4enrichment[which(l2fcol4enrichment$value >= l2fc_thres),]
    
    if (dim(filpos_log2fc_long)[1] > 0) {
      out.mac.pos <- gsub(".txt", "-enrich-MAC-positives.txt", log2fc_file)
      out.mac.pos <- paste0(output_dir, "/", out.mac.pos)

      mac.pos <- NULL
      tryCatch(
          mac.pos <- artms_enrichLog2fc(
            dataset = filpos_log2fc_long,
            specie = specie,
            heatmaps = TRUE,
            output_name = out.mac.pos,
            background = listOfGenes
          ), error = function(e){
            cat("\n\n------ (!! Error): Enrichment is not possible!\n")
            cat("                    gProfiler server is likely down\n")
            cat("                    Just wait until is up again\n\n")
            enrich = FALSE
          }
      )
      
      if(!is.null(mac.pos)){
        if (dim(mac.pos)[1] > 0) {
          write.table(
            mac.pos,
            out.mac.pos,
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE
          )
        }
      }
      
      cat("---+ Corum Protein Complex Enrichment Analysis\n")
      
      # CORUM
      positiveComplexEnriched <-
        .artms_enrichForComplexes(filpos_log2fc_long, backgroundNumber)
      
      if (dim(positiveComplexEnriched)[1] > 0) {
        out.mac.pos.corum <-
          gsub(".txt",
               "-enrich-MAC-positives-corum.txt",
               log2fc_file)
        out.mac.pos.corum <-
          paste0(output_dir, "/", out.mac.pos.corum)
        write.table(
          positiveComplexEnriched,
          out.mac.pos.corum,
          quote = FALSE,
          sep = "\t",
          row.names = FALSE,
          col.names = TRUE
        )
      }
      
      # And the heatmap
      if (dim(positiveComplexEnriched)[1] > 2) {
        out.mac.pos.corum.pdf <-
          gsub(".txt",
               "-enrich-MAC-positives-corum.pdf",
               log2fc_file)
        out.mac.pos.corum.pdf <-
          paste0(output_dir, "/", out.mac.pos.corum.pdf)
        # out.mac.pos.corum.pdf <- 'whatever.corum.positive.pdf'
        .artms_plotCorumEnrichment(
          positiveComplexEnriched,
          out.mac.pos.corum.pdf,
          "MAC+ Protein Complex Enrichment"
        )
      } else{
        cat("\t----(-) Not enough positive corum complexes to plot\n")
      }
    } else{
      cat("\t ------ Nothing is significant in the Positive site of things")
      mac.pos <- NULL
      positiveComplexEnriched <- NULL
    }
    
    
    # NEGATIVE log2fc
    cat("3) Enrichment of selected NEGATIVE significant changes\n")
    
    filneg_log2fc_long <-
      l2fcol4enrichment[which(l2fcol4enrichment$value <= -l2fc_thres),]
    
    if (dim(filneg_log2fc_long)[1] > 0) {
      out.mac.neg <- gsub(".txt", "-enrich-MAC-negatives.txt", log2fc_file)
      out.mac.neg <- paste0(output_dir, "/", out.mac.neg)
      
      mac.neg <- NULL
      tryCatch(
        mac.neg <- artms_enrichLog2fc(
          dataset = filneg_log2fc_long,
          output_name = out.mac.neg,
          specie = specie,
          heatmaps = TRUE,
          background = listOfGenes), 
        error = function(e){
          cat("\n\n------ (!! Error): Enrichment is not possible!\n")
          cat("                    gProfiler server is likely down\n")
          cat("                    Just wait until is up again\n\n")
        })
      
      if(!is.null(mac.neg)){
        if (dim(mac.neg)[1] > 0) {
          write.table(
            mac.neg,
            out.mac.neg,
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE
          )
        }
      }
      
      cat("---+ Corum Protein Complex Enrichment Analysis\n")
      
      negativesComplexEnriched <-
        .artms_enrichForComplexes(filneg_log2fc_long, backgroundNumber)
      
      if (dim(negativesComplexEnriched)[1] > 0) {
        out.mac.neg.corum <-
          gsub(".txt",
               "-enrich-MAC-negatives-corum.txt",
               log2fc_file)
        out.mac.neg.corum <-
          paste0(output_dir, "/", out.mac.neg.corum)
        write.table(
          negativesComplexEnriched,
          out.mac.neg.corum,
          quote = FALSE,
          sep = "\t",
          row.names = FALSE,
          col.names = TRUE
        )
      }
      
      # And the heatmap
      if (dim(negativesComplexEnriched)[1] > 2) {
        out.mac.neg.corum.pdf <-
          gsub(".txt",
               "-enrich-MAC-negatives-corum.pdf",
               log2fc_file)
        out.mac.neg.corum.pdf <-
          paste0(output_dir, "/", out.mac.neg.corum.pdf)
        .artms_plotCorumEnrichment(
          negativesComplexEnriched,
          out.mac.neg.corum.pdf,
          "MAC- Protein Complex Enrichment"
        )
      } else{
        cat("\t-----(-) Not enough negative corum complexes to plot\n")
      }
    } else{
      cat("\t------ Nothing is significant in the NEGATIVE site of things\n")
      mac.neg <- NULL
      negativesComplexEnriched <- NULL
    }
  } else{
    cat(">> NO ENRICHMENT of CHANGES (log2fc) SELECTED\n")
  }
  # END enrichments
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  cat(">> PRINT OUT FILES\n")
  superunified <-
    merge(abundancelongsummean,
          nbr_long,
          by = c('Bait', 'Prey'),
          all = TRUE)
  superunified <-
    merge(superunified,
          reprocondition2merge,
          by = 'Prey',
          all = TRUE)
  superunified <-
    merge(superunified,
          reprospec2merge,
          by = 'Prey',
          all = TRUE)
  
  if (grepl("ptm", isPtm)) {
    names(superunified)[grep('^Prey$', names(superunified))] <-
      'Uniprot_PTM'
    # Take the Protein ID, but being very careful about the fluomics labeling
    superunified$Prey <-
      ifelse(
        grepl("_H1N1|_H3N2|_H5N1", superunified$Uniprot_PTM),
        gsub(
          "^(\\S+?_H[1,3,5]N[1,2])_.*",
          "\\1",
          superunified$Uniprot_PTM,
          perl = TRUE
        ) ,
        gsub("^(\\S+?)_.*", "\\1", superunified$Uniprot_PTM, perl = TRUE)
      )
  }
  
  suppressMessages(superunified <-
                     artms_annotationUniprot(superunified, 'Prey', specie))
  
  # Rename (before it was just a lazy way to use another code)
  names(superunified)[grep('Bait', names(superunified))] <-
    'Condition'
  
  # ANNOTATE SPECIE
  cat("--- Annotating specie(s) in files\n")
  superunified <-
    artms_annotateSpecie(superunified, pathogen, specie)
  imputedDF <- artms_annotateSpecie(imputedDF, pathogen, specie)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # THE JITTER PLOTS
  
  if (isFluomics == TRUE) {
    cat(">> JITTERED PLOT\n")
    # Filter by number of biological replicas > 1
    superunifiedfiltered <-
      superunified[which(superunified$BioRep > 1), ]
    # Filter less than
    superunifiedfiltered <-
      superunifiedfiltered[which(superunifiedfiltered$AbMean < 30 &
                                   superunifiedfiltered$AbMean > 10),]
    # Removing carry overs
    superunifiedfiltered <-
      superunifiedfiltered[!(
        grepl("H1N1|H3N2|H5N1", superunifiedfiltered$Protein) &
          grepl("MOCK", superunifiedfiltered$Condition)
      ), ]
    
    abuJittered <-
      gsub(".txt", ".abundanceGrouped.pdf", log2fc_file)
    abuJittered <- paste0("plot.", abuJittered)
    abuJittered <- paste0(output_dir, "/", abuJittered)
    # j <- ggplot(superunifiedfiltered %>% arrange(Specie), 
    # aes(Condition,AbMean))
    # j <- j + geom_jitter(aes(colour = Specie), width = 0.3)
    if (specie == "human") {
      j <-
        ggplot(
          superunifiedfiltered %>% arrange(desc(Specie)),
          aes(
            x = Condition,
            y = AbMean,
            colour = Specie
          )
        ) #superunifiedfiltered %>% arrange(Specie)
      j <- j + geom_jitter(width = 0.3)
      j <- j + scale_colour_manual(values = c("red", "lightblue"))
    } else if (specie == "mouse") {
      j <-
        ggplot(superunifiedfiltered %>% arrange(Specie),
               aes(Condition, AbMean))
      j <- j + geom_jitter(aes(colour = Specie), width = 0.3)
      j <- j + scale_colour_manual(values = c("azure3", "red"))
    }
    j <- j + theme_minimal()
    j <-
      j + theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ))
    
    pdf(abuJittered)
    print(j)
    garbage <- dev.off()
    cat("--- done\n")
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (grepl("ptm", isPtm)) {
    cat(">> GENERATING EXTENDED DETAILED VERSION OF PH-SITE\n")
    imputedDFext <- artms_generatePhSiteExtended(
      df = imputedDF,
      pathogen = pathogen,
      specie = specie,
      ptmType = isPtm,
      output_name = paste0(output_dir, "/", log2fc_file)
    )
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cat(">> GENERATING FINAL OUTPUT FILES\n")
  if (grepl("ptm", isPtm)) {
    names(imputedDF)[grep('Protein', names(imputedDF))] <- 'Uniprot_PTM'
    imputedDF$UniprotID <- imputedDF$Uniprot_PTM
    # The virus labeling has to be taken into account 
    # when getting the uniprot id:
    imputedDF$UniprotID <-
      ifelse(
        grepl("_H1N1|_H3N2|_H5N1", imputedDF$UniprotID),
        gsub(
          "^(\\S+?_H[1,3,5]N[1,2])_.*",
          "\\1",
          imputedDF$UniprotID,
          perl = TRUE
        ) ,
        gsub("^(\\S+?)_.*", "\\1", imputedDF$UniprotID, perl = TRUE)
      )
    suppressMessages(imputedDF <-
                       artms_annotationUniprot(imputedDF, 'UniprotID', specie))
    names(imputedDF)[grep("Label", names(imputedDF))] <-
      'Comparison'
    
    imputedDF <- artms_annotateSpecie(imputedDF, pathogen, specie)
    
    # Wide version of imputedDF
    imputedDF_wide_log2fc <-
      data.table::dcast(
        data = imputedDF,
        Gene + Protein + ENTREZID + Uniprot_PTM ~ Comparison,
        value.var = 'iLog2FC',
        fill = 0
      )
    imputedDF_wide_pvalue <-
      data.table::dcast(
        data = imputedDF,
        Gene + Protein + ENTREZID + Uniprot_PTM ~ Comparison,
        value.var = 'iPvalue',
        fill = 0
      )
    
  } else if (isPtm == "global") {
    suppressMessages(imputedDF <-
                       artms_annotationUniprot(imputedDF, 'Protein', specie))
    names(imputedDF)[grep("Label", names(imputedDF))] <-
      'Comparison'
    
    # imputedDF$Specie <- ifelse(imputedDF$Protein %in% pathogen.ids$Entry, 
    # pathogen, specie)
    imputedDF <- artms_annotateSpecie(imputedDF, pathogen, specie)
    
    # Wide version of imputedDF
    imputedDF_wide_log2fc <-
      data.table::dcast(
        data = imputedDF,
        Gene + Protein + ENTREZID ~ Comparison,
        value.var = 'iLog2FC',
        fill = 0
      )
    imputedDF_wide_pvalue <-
      data.table::dcast(
        data = imputedDF,
        Gene + Protein + ENTREZID ~ Comparison,
        value.var = 'iPvalue',
        fill = 0
      )
  } else{
    stop("\nWRONG isPTM SELECTED. 
         OPTIONS AVAILABLE: global, ptmph, yesphsite\n")
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # boxplot of relative abundances
  cat(">> PLOT OUT: TOTAL NUMBER OF PROTEINS/SITES QUANTIFIED\n")
  numimputedfinal <-
    gsub(".txt", ".TotalNumberQuantifications.pdf", log2fc_file)
  numimputedfinal <- paste0("plot.", numimputedfinal)
  numimputedfinal <- paste0(output_dir, "/", numimputedfinal)
  pdf(numimputedfinal)
  .artms_plotNumberProteinsImputedLog2fc(imputedDF)
  garbage <- dev.off()
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (isPtm == "global") {
    cat(">> CLUSTERING ANALYSIS OF QUANTIFICATIONS\n")
    
    # GET THE LIST OF SIGNIFICANTS FOR THE EXPERIMENT(S)
    list_of_significants <-
      unique(imputedDF$Protein[which(abs(imputedDF$iLog2FC > 1) &
                                       imputedDF$iPvalue < 0.05)])
    
    # AND APPLY THE FILTER
    data.select <-
      imputedDF[which(imputedDF$Protein %in% list_of_significants), ]
    
    # GENE BASED -> heatmaps
    hasdc <-
      data.table::dcast(
        data = data.select[which(data.select$imputed == "no"), ],
        Gene + Protein ~ Comparison,
        value.var = "iLog2FC",
        fun.aggregate = median,
        fill = 0
      )
    
    # EXPERIMENT BASED -> PCA
    hasdcexp <-
      data.table::dcast(
        data = data.select[which(data.select$imputed == "no"), ],
        Comparison ~ Gene + Protein,
        value.var = "iLog2FC",
        fun.aggregate = median,
        fill = 0
      )
    
    # ----------------------------------------------------------------------
    # CLUSTERING ANALYSIS
    
    # GENE BASED
    rownames(hasdc) <- paste0(hasdc$Gene, "_", hasdc$Protein)
    vamos <- within(hasdc, rm(Gene, Protein))
    
    # if this dataset only have one comparison,
    # this analysis does not makes sense: check it out:
    if (dim(vamos)[2] > 1) {
      venga <- as.matrix(vamos)
      
      # EXPERIMENT BASED
      rownames(hasdcexp) <- hasdcexp$Comparison
      vamosexp <- within(hasdcexp, rm(Comparison))
      vengaexp <- as.matrix(vamosexp)
      
      # PCA AND CORRELATION ANALYSIS
      cat("--- Correlation plots\n")
      df.cor.matrix <-
        round(cor(venga, use = "pairwise.complete.obs"), 2)
      file_corr_l2fc <- gsub(".txt", ".log2fc-corr.pdf", log2fc_file)
      file_corr_l2fc <- paste0(output_dir, "/", file_corr_l2fc)
      pdf(file_corr_l2fc, width = 12, height = 9)
      corrplot::corrplot(
        df.cor.matrix,
        type = "upper",
        tl.pos = "td",
        method = "circle",
        tl.cex = 0.9,
        tl.col = 'black',
        tl.srt = 45,
        # order = "hclust",
        diag = TRUE
      )
      PerformanceAnalytics::chart.Correlation(venga,
                                    histogram = TRUE,
                                    pch = 25,
                                    main = "Correlation between Comparisons")
      garbage <- dev.off()
      
      # BASED ON GROUPS
      pca.hasdcexp <-
        FactoMineR::PCA(
          hasdcexp[, -c(1)],
          scale.unit = FALSE,
          ncp = 4,
          graph = FALSE
        )
      
      pca_all <- factoextra::fviz_pca_ind(
        pca.hasdcexp,
        labelsize = 3,
        repel = TRUE,
        habillage = as.factor(hasdcexp$Comparison),
        addEllipses = FALSE,
        ellipse.level = 0.95
      )
      
      cat("--- PCA, individuals plot\n")
      file_pca_l2fc <-
        gsub(".txt", ".log2fc-individuals-pca.pdf", log2fc_file)
      file_pca_l2fc <- paste0(output_dir, "/", file_pca_l2fc)
      pdf(file_pca_l2fc, width = 9, height = 7)
      print(pca_all)
      garbage <- dev.off()
      
      # Determine the OPTIMAL NUMBER OF CLUSTERS:
      
      # Elbow method
      e1 <-
        factoextra::fviz_nbclust(venga, kmeans, method = "wss") +
        geom_vline(xintercept = 4, linetype = 2) +
        labs(subtitle = "kmeans Elbow method")
      e2 <-
        factoextra::fviz_nbclust(venga, cluster::pam, method = "wss") +
        geom_vline(xintercept = 4, linetype = 2) +
        labs(subtitle = "PAM Elbow method")
      
      # Silhouette method
      k1 <-
        factoextra::fviz_nbclust(venga, kmeans, method = "silhouette") +
        labs(subtitle = "kmeans Silhouette method")
      k2 <-
        factoextra::fviz_nbclust(venga, cluster::pam, method = "silhouette") +
        labs(subtitle = "pam Silhouette method")
      
      
      # Create a dendrogram
      cat("--- Dendrogram\n")
      res.dist <-
        factoextra::get_dist(vamosexp, stand = TRUE, method = "minkowski")
      hc <- hclust(res.dist)
      file_dendro_l2fc <-
        gsub(".txt", ".log2fc-dendro.pdf", log2fc_file)
      file_dendro_l2fc <- paste0(output_dir, "/", file_dendro_l2fc)
      pdf(file_dendro_l2fc, width = 9, height = 7)
      plot(hc)
      garbage <- dev.off()
      
      # COMPLEXHEATMAP Heatmap with a specified number of optimal clusters
      n = 10
      pam.res <- pam(vamos, k = n)
      
      cp1 <- factoextra::fviz_cluster(pam.res)
      cp2 <-
        factoextra::fviz_silhouette(pam.res, print.summary = FALSE)
      
      cat("--- Plots to determine optimal number of clusters\n")
      file_clusterplots_l2fc <-
        gsub(".txt", ".log2fc-clusters.pdf", log2fc_file)
      file_clusterplots_l2fc <-
        paste0(output_dir, "/", file_clusterplots_l2fc)
      pdf(file_clusterplots_l2fc,
          width = 9,
          height = 7)
      print(e1)
      print(e2)
      print(k1)
      print(k2)
      print(cp1)
      print(cp2)
      garbage <- dev.off()
      
      cat("--- Cluster heatmaps (10 clusters)\n")
      hmap <- ComplexHeatmap::Heatmap(
        vamos,
        name = paste0("Clusters ", "(n = ", n, ")"),
        col = circlize::colorRamp2(c(-3, 0, 3), c(
          "firebrick1", "black", "olivedrab1"
        )),
        heatmap_legend_param = list(
          color_bar = "continuous",
          legend_direction = "horizontal",
          legend_width = unit(5, "cm"),
          title_position = "topcenter",
          title_gp = gpar(fontsize = 15, fontface = "bold")
        ),
        split = paste0("", pam.res$clustering),
        row_title = "Genes",
        row_title_side = "left",
        row_title_gp = gpar(fontsize = 15, fontface = "bold"),
        show_row_names = FALSE,
        column_title = "Relative Quantifications",
        column_title_side = "top",
        column_title_gp = gpar(fontsize = 10, fontface = "bold"),
        column_title_rot = 0,
        show_column_names = TRUE,
        cluster_columns = FALSE,
        clustering_distance_columns = function(x)
          as.dist(1 - cor(t(x))),
        clustering_method_columns = "ward.D2",
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "ward.D2",
        row_dend_width = unit(30, "mm"),
        column_dend_height = unit(30, "mm"),
        # top_annotation=colAnn,
        top_annotation_height = unit(1.75, "cm"),
        # bottom_annotation=sampleBoxplot,
        bottom_annotation_height = unit(4, "cm"),
        column_names_gp = gpar(fontsize = 10)
      )
      file_clusterheat_l2fc <-
        gsub(".txt", ".log2fc-clusterheatmap.pdf", log2fc_file)
      file_clusterheat_l2fc <-
        paste0(output_dir, "/", file_clusterheat_l2fc)
      pdf(file_clusterheat_l2fc,
          width = 12,
          height = 10)
      ComplexHeatmap::draw(hmap,
                           heatmap_legend_side = "top",
                           annotation_legend_side = "right")
      garbage <- dev.off()
      
      cat("--- Enrichment analysis of the clusters\n")
      cl_number <- pam.res$clustering
      dfclusters <- as.data.frame(cl_number)
      dfclusters$ids <- row.names(dfclusters)
      dfclusters$Gene <- gsub("(.*)(_)(.*)", "\\1", dfclusters$ids)
      dfclusters$Protein <- gsub("(.*)(_)(.*)", "\\3", dfclusters$ids)
      
      # Making sure we have unique genes in each comparison 
      # (the PTM might bring redundancy)
      pretmp <- dfclusters[c('Gene', 'cl_number')]
      pretmp <- unique(pretmp)
      
      tmp <- split(pretmp$Gene, pretmp$cl_number, drop = TRUE)
      
      if (specie == "human") {
        enrichgenes <-
          artms_enrichProfiler(
            tmp,
            categorySource = c(
              'GO:BP',
              'GO:MF',
              'GO:CC',
              'KEGG',
              'REAC',
              'CORUM',
              'HPA',
              'OMIM'
            ),
            specie = 'hsapiens',
            listOfGenes
          ) # 'HP'
      } else if (specie == "mouse") {
        enrichgenes <-
          artms_enrichProfiler(
            tmp,
            categorySource = c('GO:BP', 
                               'GO:MF', 
                               'GO:CC', 
                               'KEGG', 
                               'REAC', 
                               'CORUM'),
            specie = 'mmusculus',
            listOfGenes
          )
      } else{
        stop(
          "\n\n\ntOhhh no, this specie",
          specie,
          " is not supported in the enrichment!!\n\n\n"
        )
      }
      
      file_clusterheatenrich_l2fc <-
        gsub(".txt",
             ".log2fc-clusterheatmap-enriched.txt",
             log2fc_file)
      file_clusterheatenrich_l2fc <-
        paste0(output_dir, "/", file_clusterheatenrich_l2fc)
      write.table(
        enrichgenes,
        file_clusterheatenrich_l2fc,
        col.names = TRUE,
        row.names = FALSE,
        sep = "\t",
        quote = FALSE
      )
      
      file_clusterheatdata_l2fc <-
        gsub(".txt", ".log2fc-clusterheatmap.txt", log2fc_file)
      file_clusterheatdata_l2fc <-
        paste0(output_dir, "/", file_clusterheatdata_l2fc)
      write.table(
        dfclusters,
        file_clusterheatdata_l2fc,
        col.names = TRUE,
        row.names = FALSE,
        sep = "\t",
        quote = FALSE
      )
    }
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cat(">> WRITTING ALL THE OUTPUT FILES\n")
  
  # PRINT OUT IMPUTED
  outlog2fcImpute <- gsub(".txt", "-log2fc-long.txt", log2fc_file)
  outlog2fcImpute <- paste0(output_dir, "/", outlog2fcImpute)
  write.table(
    imputedDF,
    outlog2fcImpute,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  
  outlog2fc <- gsub(".txt", "-log2fc-wide.txt", log2fc_file)
  outlog2fc <- paste0(output_dir, "/", outlog2fc)
  write.table(
    log2fc_file_splc,
    outlog2fc,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  
  outmodeqcLong <- gsub(".txt", "-abundance-long.txt", log2fc_file)
  outmodeqcLong <- paste0(output_dir, "/", outmodeqcLong)
  write.table(
    superunified,
    outmodeqcLong,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  
  outmodelqc <- gsub(".txt", "-abundance-wide.txt", log2fc_file)
  outmodelqc <- paste0(output_dir, "/", outmodelqc)
  write.table(
    modelqc_file_splc,
    outmodelqc,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  
  outexcel <- gsub(".txt", "-summary.xlsx", log2fc_file)
  outexcel <- paste0(output_dir, "/", outexcel)
  
  if (enrich) {
    # But now check whether is a PTM case:
    if (grepl("ptm", isPtm)) {
      list_of_datasets <- list(
        "log2fcImputed" = imputedDF,
        "log2fcImpExt" = imputedDFext,
        "wide_iLog2fc" = imputedDF_wide_log2fc,
        "wide_iPvalue" = imputedDF_wide_pvalue,
        "enrichALL" = mac.allsig,
        "enrichMACpos" = mac.pos,
        "enrichMACneg" = mac.neg,
        "enMACallCorum" = allsigComplexEnriched,
        "enMACposCorum" = positiveComplexEnriched,
        "enMACnegCorum" = negativesComplexEnriched
      )
    } else if (grepl("global", isPtm)) {
      list_of_datasets <- list(
        "log2fcImputed" = imputedDF,
        "wide_iLog2fc" = imputedDF_wide_log2fc,
        "wide_iPvalue" = imputedDF_wide_pvalue,
        "enrichALL" = mac.allsig,
        "enrich-MACpos" = mac.pos,
        "enrich-MACneg" = mac.neg,
        "enMACallCorum" = allsigComplexEnriched,
        "enMACposCorum" = positiveComplexEnriched,
        "enMACnegCorum" = negativesComplexEnriched
      )
      
    } else{
      stop("Oh no!! This will fail if you are using UB!!\n")
    }
  } else if (!enrich) {
    cat("-----(-) Enrichment was not selected\n")
    if (grepl("ptm", isPtm)) {
      list_of_datasets <- list(
        "log2fcImputed" = imputedDF,
        "log2fcImpExt" = imputedDFext,
        "wide_iLog2fc" = imputedDF_wide_log2fc,
        "wide_iPvalue" = imputedDF_wide_pvalue
      )
    } else if (grepl("global", isPtm)) {
      list_of_datasets <- list(
        "log2fcImputed" = imputedDF,
        "wide_iLog2fc" = imputedDF_wide_log2fc,
        "wide_iPvalue" = imputedDF_wide_pvalue
      )
    } else{
      stop("Oh no!! This will fail if you are using UB!!\n")
    }
  } else{
    stop(
      "\n\nYOU SHOULD NEVER SEE THIS MESSAGE. 
      IF you do, dude, check the source code urgently\n\n"
    )
    # The script should have crashed by this point. 
    # If it gets up to here... it would be very weird
  }
  
  # Defining style for the header
  hs <-
    openxlsx::createStyle(
      fontName = "Arial",
      fontColour = "white",
      fgFill = "#000000",
      textDecoration = "Bold",
      border = "Bottom"
    )
  openxlsx::write.xlsx(
    list_of_datasets,
    file = outexcel,
    asTable = TRUE,
    headerStyle = hs
  )
  
  cat("--- OUTPUT FILES IN FOLDER <", output_dir, ">\n")
  cat("- EXCEL: ", outexcel, "\n")
  cat("- Log2fc Wide: ", outlog2fc, "\n")
  cat("- Log2fc Impute: ", outlog2fc, "\n")
  cat("- AbundanceLong: ", outmodeqcLong, "\n")
  cat("- AbundanceWide: ", outmodelqc, "\n")
  
  if (enrich == TRUE) {
    cat("- ENRICHMENT files should also be out\n")
  }
  
  cat("\n>> SUPER ANALYSIS IS DONE\n\n")
}


# ------------------------------------------------------------------------------
#' @title Adding a column with the specie name
#'
#' @description Adding the specie name to every protein.
#' This makes more sense if there are more than one specie in the dataset,
#' which must be specified in the `pathogen` option. Influenza is a special
#' case that it does not need to be specified, as far as the proteins were
#' originally annotated as `INFLUENZAGENE_STRAIN`
#' (strains covered `H1N1`, `H3N2`, `H5N1`), as for example, `NS1_H1N1`
#' @param df (data.frame) with a `Protein` column (of uniprot ids)
#' @param pathogen (char) Is there a pathogen in the dataset as well?
#' if it does not, then use `pathogen = nopathogen` (default). Supported`tb` 
#' (Tuberculosis),
#' `lpn` (Legionella)
#' @param specie (char) Host organism (supported for now: `human` or `mouse`)
#' @return (data.frame) The same data.frame but with an extra column 
#' specifying the specie
#' @keywords annotation, specie
#' @examples
#' # Adding a new column with the main specie of the data. Easy.
#' # But the main functionality is to add both the host-specie and a pathogen,
#' # which is not illustrated in this example
#' data_with_specie <- artms_annotateSpecie(df = artms_data_ph_msstats_results,
#'                                          specie = "human")
#' @export
artms_annotateSpecie <- function(df,
                                 pathogen = "nopathogen",
                                 specie) {
  if (pathogen == "nopathogen") {
    # Influenza is treated differently
    df$Specie <-
      ifelse(grepl("_H1N1|_H3N2|_H5N1", df$Protein),
             "Influenza",
             specie)
  } else if (pathogen == "tb") {
    pathogen.ids <- artms_data_pathogen_TB
    df$Specie <-
      ifelse(df$Protein %in% pathogen.ids$Entry, pathogen, specie)
  } else if (pathogen == "lpn") {
    pathogen.ids <- artms_data_pathogen_LPN
    df$Specie <-
      ifelse(df$Protein %in% pathogen.ids$Entry, pathogen, specie)
  }
  return(df)
}

# ------------------------------------------------------------------------------
#' @title Generate ph-site specific detailed file
#'
#' @description Generate extended detailed ph-site file, where every line is a
#' ph site instead of a peptide. Therefore, if one peptide has multiple ph sites
#' it will be breaking down in each of the sites. This file will help generate
#' input files for tools as [Phosfate](http://phosfate.com/) or
#' [PHOTON](https://github.com/jdrudolph/photon)
#' @param df (data.frame) of log2fc and imputed values
#' @param pathogen (char) Is there a pathogen in the dataset as well? Available
#' pathogens are `tb` (Tuberculosis), `lpn` (Legionella). If it is not,
#' then use `nopathogen` (default).
#' @param specie (char) Main organism (supported for now: `human` or `mouse`)
#' @param ptmType (char) It must be a ptm-site quantification dataset. Either:
#' yes: `ptmsites` (for site specific analysis), or
#' `ptmph` (Jeff's script output evidence file).
#' @param output_name (char) A output file name (extension `.txt` required)
#' @return (data.frame) extended version of the ph-site
#' @keywords external, tools, phosfate
#' @examples \donttest{
#' artms_generatePhSiteExtended(df = dfobject, 
#'                              specie = "mouse", 
#'                              ptmType = "ptmsites",
#'                              output_name = log2fc_file)
#' }
#' @export
artms_generatePhSiteExtended <-
  function(df, 
           pathogen = "nopathogen", 
           specie, 
           ptmType,
           output_name) {
    
    imputedDFext <- NULL
    
    if (ptmType == "ptmph") {
      imputedDFext <- df
      names(imputedDFext)[grep('^Protein$', names(imputedDFext))] <-
        'Uniprot_PTM'
      # Take the Protein ID, but being very careful about the fluomics labeling
      imputedDFext$Protein <- ifelse(
        grepl("_H1N1|_H3N2|_H5N1", imputedDFext$Uniprot_PTM),
        gsub(
          "^(\\S+?_H[1,3,5]N[1,2])_.*",
          "\\1",
          imputedDFext$Uniprot_PTM,
          perl = TRUE
        ) ,
        gsub("^(\\S+?)_.*", "\\1", imputedDFext$Uniprot_PTM, perl = TRUE)
      )
      
      # Extract sites from Uniprot_PTM
      imputedDFext$PTMsite <-
        gsub("^(\\S+?)(_ph.*)", "\\2", imputedDFext$Uniprot_PTM, perl = TRUE)
      imputedDFext$PTMsite <- gsub("^(_ph)", "", imputedDFext$PTMsite)
      imputedDFext$PTMsite <- gsub("_ph", ",", imputedDFext$PTMsite)
      # And create independent columns for each of them
      imputedDFext <-
        imputedDFext %>% mutate(PTMsite = strsplit(PTMsite, ",")) %>% tidyr::unnest(PTMsite)
      
      suppressMessages(imputedDFext <-
                         artms_annotationUniprot(imputedDFext, 
                                                 'Protein', 
                                                 specie))
      names(imputedDFext)[grep("^Label$", names(imputedDFext))] <-
        'Comparison'
      
      imputedDFext <-
        artms_annotateSpecie(imputedDFext, pathogen, specie)
    } else if (ptmType == "ptmsites") {
      imputedDFext <- df
      #1. Change the Protein name
      names(imputedDFext)[grep('^Protein$', names(imputedDFext))] <- 
        'Uniprot_PTM'
      
      # 2. Make a copy of Uniprot_PTM to operate on it
      imputedDFext$PTMone <- imputedDFext$Uniprot_PTM
      
    # 3. Create independent columns for each of them
    imputedDFext <-
      imputedDFext %>% mutate(PTMone = strsplit(PTMone, ",")) %>% tidyr::unnest(PTMone)
      
      # 4. And take the labels:
      imputedDFext$Protein <-
        ifelse(
          grepl("_H1N1|_H3N2|_H5N1", imputedDFext$PTMone),
          gsub("^(\\S+?_H[1,3,5]N[1,2])_.*", "\\1",
               imputedDFext$PTMone, perl = TRUE
          ) ,
          gsub("^(\\S+?)_.*", "\\1", imputedDFext$PTMone, perl = TRUE)
        )
      imputedDFext$PTMAA <-
        gsub("(\\S+)(_[S,T,Y,K])(\\d+)", "\\2", imputedDFext$PTMone)
      imputedDFext$PTMsite <-
        gsub("(\\S+)(_[S,T,Y,K])(\\d+)", "\\3", imputedDFext$PTMone)
      imputedDFext <-
        artms_annotationUniprot(imputedDFext, 'Protein', specie)
      names(imputedDFext)[grep("^Label$", names(imputedDFext))] <-
        'Comparison'
      
      # imputedDFext$Specie <- ifelse(grepl("_H1N1|_H3N2|_H5N1", 
      # imputedDFext$Protein), "Influenza", specie)
      # imputedDFext$Specie <- ifelse(imputedDFext$Protein 
      # %in% pathogen.ids$Entry, pathogen, specie)
      imputedDFext <-
        artms_annotateSpecie(imputedDFext, pathogen, specie)
    } else{
      stop(
        "--- (!!!) Only 'ptmph' or 'ptmsites' allowed for argument <ptmType>\n"
      )
    }
    outlog2fcImputext <-
      gsub(".txt", "-imputedL2fcExtended.txt", output_name)
    write.table(
      imputedDFext,
      outlog2fcImputext,
      quote = FALSE,
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE
    )
    return(imputedDFext)
    cat("--- ph-site extended version ready\n")
  }


# ------------------------------------------------------------------------------
# @title Imputing missing values
#
# @description When a value is completely missed in one of the conditions,
# the `log2fc = Inf / -Inf`. This function imputes those values, i.e.,
# will assign 'artificial' values.
# The imputation method works as follow. The assumption is that those
# proteins are likely present as well in those conditions where are found
# missed, but due to the small sampling (usually 2 or 3 biological replicas)
# and other proteomics related issues, those proteins didn't make it through
# the level of detection.
# Therefore, a small intensity (sampled from the bottom 10 intensity values)
# will be assigned to the protein/site in the missing condition,
# and the new log2fc is re-calculated out of the MSstats box.
# Two issues are addressed as follows:
# 1. If a protein has been consistently identified in one of the conditions,
# it will stay
# 2. But if the intensity value in those conditions was too low,
# then the log2fc will be also low
# @param dflog2fcinfinites data.frame of proteins with only infinite
# values from the msstats results file
# @param dfmq Abundance data, which will be used to know the details of
# reproducibility
# @return Imputed missing values
# @keywords internal, imputation, log2fc, quantifications, missing values
.artms_imputeMissingValues <- function(dflog2fcinfinites, dfmq) {
  # The comparsions
  contrast <- unique(dflog2fcinfinites$Label)
  
  # Select the IDs to impute
  ids2impute <- unique(dflog2fcinfinites$Protein)
  
  # Take the abundance values for all the proteins
  abu2imp <- .artms_loadModelqcBasic(dfmq)
  # Aggregate the technical replica by choosing the maximum value
  abu2imp2 <-
    aggregate(Abundance ~ Protein + Condition + BioReplicate,
              data = abu2imp,
              FUN = mean)
  
  # Check things that will be imputed
  # dfdc.ni <- data.table::dcast(data=abu2imp2, 
  # Protein~BioReplicate, value.var = "Abundance")
  
  # Two possible options here.
  # 1. Select based on the bottom x%
  # # Imputing the missing values by selecting randomly from the bottom 5%
  # theMin <- min(dfmq$ABUNDANCE)
  # # Select the 5% quartile as the maximum value to sample from
  # theMax <- quantile(dfmq$ABUNDANCE, probs = .05)
  #
  # 2. Select the bottom 20 intensities
  # Grab the bottom 30 intensities in the dataset
  dfmqOrdered <- dfmq[order(dfmq$ABUNDANCE, decreasing = FALSE), ]
  
  numberFromBottom <- 10
  abuBottom <- head(dfmqOrdered$ABUNDANCE, n = numberFromBottom)
  
  theMin <- abuBottom[1]
  theMax <- abuBottom[numberFromBottom]
  
  # Generating the numbers from which we are going to sample
  numbers2sample <- seq(from = theMin, to = theMax, by = .00001)
  
  # dcast on abundance and fill with random numbers between the minimum and q05
  suppressWarnings(
    dfdc.im <-
      data.table::dcast(
        data = abu2imp2,
        Protein ~ BioReplicate,
        value.var = "Abundance",
        fill = sample(numbers2sample, replace = FALSE)
      )
  )
  
  # Needs to aggregate on biological replicas
  # 1. Melt on biological replicas
  dfdc.melt <-
    reshape2::melt(
      dfdc.im,
      id.vars = c('Protein'),
      value.name = 'Abundance',
      variable.name = 'BioReplicate'
    )
  # 2. Get the condition
  dfdc.melt$Condition <-
    gsub("(.*)(-)(.*)", "\\1", dfdc.melt$BioReplicate)
  # 3. Dcast and Aggregate on the condition, taking the mean
  dfdc.final <-
    data.table::dcast(
      data = dfdc.melt,
      Protein ~ Condition,
      value.var = 'Abundance',
      fun.aggregate = mean
    )
  # 4. Filter by proteins to impute
  dfdc.final <-
    dfdc.final[which(dfdc.final$Protein %in% ids2impute), ]
  
  for (c in contrast) {
    # cat("\t",c," --> ")
    x <- gsub("(.*)(-)(.*)", "\\1", c)
    y <- gsub("(.*)(-)(.*)", "\\3", c)
    # cat("log2fc(",x, " - ", y,")\n")
    
    # Renaming the comparision name just for illustration purposes
    rnc <- paste0("l2fc_", c)
    
    dfdc.final[[rnc]] <- dfdc.final[[x]] - dfdc.final[[y]]
  }
  
  # Select only log2fc columns
  imputedL2FCValues <-
    dfdc.final[grepl("Protein|l2fc_", colnames(dfdc.final))]
  
  # Melt again
  imputedL2FCmelted <-
    reshape2::melt(
      imputedL2FCValues,
      id.vars = c('Protein'),
      variable.name = 'Label',
      value.name = 'iLog2FC'
    )
  # Now let's get it ready for merging with the values to be 
  # imputed at dflog2fcinfinites
  imputedL2FCmelted$Label <-
    gsub("l2fc_", "", imputedL2FCmelted$Label)
  
  # And let's add p-values
  samplingPvalue <- seq(from = 0.01, to = 0.05, by = .0000001)
  # And add imputed pvalues
  imputedL2FCmelted$iPvalue <-
    sample(samplingPvalue,
           size = nrow(imputedL2FCmelted),
           replace = FALSE)
  
  return (imputedL2FCmelted)
}

# ------------------------------------------------------------------------------
# @title Load limited columns from abundance (modelqc) annotated
#
# @description Load limited columns from abundance (modelqc) annotated
# @param df_input (data.frame) with the raw abundance data (modelqc)
# @param specie (char) Specie name for annotation purposes
# @param ptmis (char) Specify whether is a PTM dataset: `global`, `ptmsites`,
# `ptmph`
# @return annotated data.frame of abundance data
# @keywords abundance, annotated
.artms_loadModelQCstrict <- function (df_input, specie, ptmis) {
  cat("--- Loading abundance for proteins found in all biological replicas\n")
  cat("--- With respect to the ptm: ", ptmis)
  # Remove empty entries
  if (any(df_input$PROTEIN == "")) {
    df_input <- df_input[-which(df_input$PROTEIN == ""), ]
  }
  df_input$PROTEIN <-
    gsub("(^sp\\|)(.*)(\\|.*)", "\\2", df_input$PROTEIN)
  df_input$PROTEIN <- gsub("(.*)(\\|.*)", "\\1", df_input$PROTEIN)
  
  # Technical replicas: aggregate on the mean the technical replicas
  b <-
    aggregate(
      ABUNDANCE ~ PROTEIN + GROUP_ORIGINAL + SUBJECT_ORIGINAL,
      data = df_input,
      FUN = mean
    )
  
  datadc <-
    data.table::dcast(
      data = b,
      PROTEIN ~ GROUP_ORIGINAL,
      value.var = 'ABUNDANCE',
      fun.aggregate = mean
    )
  
  names(datadc)[grep('PROTEIN', names(datadc))] <- 'Protein'
  if (grepl("ptm", ptmis)) {
    # if is a PTM dataset we don't need the real gene names for now,
    # we need to use the Uniprot_ptm notation
    datadc$Gene <- datadc$Protein
    send_back <- datadc
  } else{
    suppressMessages(send_back <-
                       artms_annotationUniprot(datadc, 'Protein', specie))
  }
  return(send_back)
}


#------------------------------------------------------------------------------
# @title Load the basic ModelQC file
#
# @param data (data.frame) of the ModelQC file
# @return (data.frame) of the modelqc file with the columns Protein, Abundance,
# Condition, BioReplicate
# @keywords internal, loading
.artms_loadModelqcBasic <- function(data) {
  if (length(grep(";", data$PROTEIN)) > 0)
    data <-
      data[-grep(";", data$PROTEIN), ]
  if ("PROTEIN" %in% colnames(data)) {
    names(data)[grep("PROTEIN", names(data))] <- 'Protein'
  } else{
    cat("ERROR: you should check the abundance file
        because something is seriously wrong!\n")
    stop("Abort mission\n!")
  }
  if ("ABUNDANCE" %in% colnames(data)) {
    names(data)[grep("ABUNDANCE", names(data))] <- 'Abundance'
  } else{
    cat("ERROR: you should check the abundance file 
        because something is seriously wrong!\n")
    stop("Abort mission\n!")
  }
  if ("GROUP_ORIGINAL" %in% colnames(data)) {
    names(data)[grep("GROUP_ORIGINAL", names(data))] <- 'Condition'
  } else{
    cat("ERROR: you should check the abundance file 
        because something is seriously wrong!\n")
    stop("Abort mission\n!")
  }
  if ("SUBJECT_ORIGINAL" %in% colnames(data)) {
    names(data)[grep("SUBJECT_ORIGINAL", names(data))] <- 'BioReplicate'
  } else{
    cat("ERROR: you should check the abundance file 
        because something is seriously wrong!\n")
    stop("Abort mission\n!")
  }
  data <-
    subset(data, select = c(Protein, Abundance, Condition, BioReplicate))
  return(data)
}

# ------------------------------------------------------------------------------
# @title Merge abundance and number of biological replicates per condition
#
# @description Merge abundance and number of biological replicates
# per condition
# @param df_input (data.frame) Abundance input file
# @param repro (data.frame) Reproducibility data.frame
# @param specie (char) Specie for annotation purposes
# @return (data.frame) of abundance merged with reproducibility info
# @keywords abundance, reproducibility, merging
.artms_mergeAbNbr <- function (df_input, repro, specie) {
  # Remove empty entries
  if (any(df_input$PROTEIN == "")) {
    df_input <- df_input[-which(df_input$PROTEIN == ""), ]
  }
  df_input$PROTEIN <-
    gsub("(^sp\\|)(.*)(\\|.*)", "\\2", df_input$PROTEIN)
  df_input$PROTEIN <- gsub("(.*)(\\|.*)", "\\1", df_input$PROTEIN)
  
  # TECHNICAL REPLICAS: if there are technical replicas, 
  # this means that we will find
  # two values for the same protein in the same bioreplica, 
  # therefore we need to aggregate first just in case:
  df_input <-
    aggregate(
      ABUNDANCE ~ PROTEIN + GROUP_ORIGINAL + SUBJECT_ORIGINAL,
      data = df_input,
      FUN = mean
    )
  
  dc_input <-
    data.table::dcast(
      data = df_input[, c('PROTEIN', 'ABUNDANCE', 'GROUP_ORIGINAL')],
      PROTEIN ~ GROUP_ORIGINAL,
      value.var = 'ABUNDANCE',
      fun.aggregate = mean,
      fill = 0
    )
  names(dc_input)[grep('PROTEIN', names(dc_input))] <- 'Protein'
  
  colnames(repro) <- paste("NumBR", colnames(repro), sep = "_")
  colnames(repro)[1] <- 'Protein'
  dc_input <- merge(dc_input, repro, by = c('Protein'))
  
  return(dc_input)
}

# ------------------------------------------------------------------------------
# @title Merge changes (log2fc) and number of biological replicates per
# condition
#
# @description Merge changes, i.e., MSstats results file of quantified changes,
# with the number of biological replicates per condition
# @param df_input Changes data.frame
# @param repro Reproducibility data.frame
# @param specie Specie for annotation purposes
# @return Merged data.frame of changes and reproducibility information
# @keywords changes, log2fc, reproducibility, merging
.artms_mergeChangesNbr <- function (df_input, repro, specie) {
  # # Remove the weird empty proteins
  # if(any(df_input$Protein == "")){ df_input <- 
  # df_input[-which(df_input$Protein == ""),]}
  # df_input$Protein <- gsub("(^sp\\|)(.*)(\\|.*)", "\\2", df_input$Protein )
  # df_input$Protein <- gsub("(.*)(\\|.*)", "\\1", df_input$Protein )
  
  input_melt = reshape2::melt(data = df_input[, c('Protein', 
                                                  'Label', 
                                                  'log2FC', 
                                                  'adj.pvalue'), ], 
                              id.vars = c('Protein', 'Label'))
  input_dcast = data.table::dcast(Protein ~ Label + variable,
                                  data = input_melt,
                                  value.var = c('value'))
  
  colnames(repro) <- paste("NumBR", colnames(repro), sep = "_")
  colnames(repro)[1] <- 'Protein'
  input_dcast <- merge(input_dcast, repro, by = c('Protein'))
  
  # Move Gene name to the left:
  return(input_dcast)
}

# ------------------------------------------------------------------------------
# @title Plot the total number of quantified proteins in each condition
#
# @description
# @keys internal, plot, counts
# @param (data.frame) Data frame of imputed log2fc
.artms_plotNumberProteinsImputedLog2fc <- function(data) {
  x <- data[c('Protein', 'Comparison')]
  y <- unique(x)
  z <- ggplot(y, aes(x = Comparison, fill = Comparison))
  z <- z + geom_bar(stat = "count")
  z <-
    z + theme(axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ))
  z <-
    z + geom_text(
      stat = 'count',
      aes(label = ..count..),
      vjust = -0.5,
      size = 2.7
    )
  z <- z + ggtitle("Unique Proteins in Comparisons")
  print(z)
}

# ------------------------------------------------------------------------------
# @title Filter: Remove proteins below some threshold of minimal reproducibility
#
# @description If a protein is not found in a minimal number of
# biological replicates in at least one of the conditions, it is removed
# @param dfi (data.frame) Data.frame with biological replicates information
# @param mnbr (int) minimal number of biological replicates
# @return (data.frame) a filtered `dfi`
# @keywords internal, filter, bioreplicates, reproducibility
.artms_RemoveProtBelowThres <- function(dfi, mnbr) {
  theComparisons2check <- unique(dfi$Label)
  for (onlyonecomp in (theComparisons2check)) {
    ax <- gsub("(.*)(-)(.*)", "\\1", onlyonecomp)
    ay <- gsub("(.*)(-)(.*)", "\\3", onlyonecomp)
    
    # If the condition is not met, i.e., if all the proteins are
    # found in at least X biological replicas, then
    # it would remove the whole thing.
    if (dim(dfi[which((dfi[[ax]] < mnbr) &
                      (dfi[[ay]] < mnbr)), ])[1] > 0) {
      dfi <-
        dfi[-which(((dfi[[ax]] < mnbr) &
                      (dfi[[ay]] < mnbr)) & (dfi$Label == onlyonecomp)), ]
    }
  }
  return(dfi)
}

# ------------------------------------------------------------------------------
# @title Select and label the condition more abundant in a quantification
#
# @description Select and label the condition more abundant in a quantification
# - If log2fc > 0 the condition on the left ('numerator') is the most abundant
# - If log2fc < 0 the condition on the right ('denominator') 
# is the most abundant
# @param a (char) log2fc column
# @param b (char) comparison column
# @return (char) One of the conditions from the comparison
# @keywords internal, selection, labeling
.artms_selectTheOneLog2fc <- function(a, b) {
  thrs <- 0
  # a should be the log2fc column
  if (a > thrs) {
    sb <- gsub("(.*)(-)(.*)", "\\1", b)
  } else if (a < -thrs) {
    sb <- gsub("(.*)(-)(.*)", "\\3", b)
  } else {
    sb <- 'NA'
  }
  return(sb)
}
