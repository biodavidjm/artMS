# CREATE SAMPLE DATA FILES FOR THE PACKAGE----

setwd('~/github/biodavidjm/artMS/')


# CONFIGURATION FILE----

# library(yaml)
artms_config <- yaml.load_file("~/github/biodavidjm/artMS/inst/extdata/artms_config.yaml")
save(artms_config, file = '~/github/biodavidjm/artMS/data/artms_config.RData', compress = TRUE)

# GENERATE RANDOM FILE----
artms_data_randomDF <- data.frame(replicate(10, sample(0:1, 100, rep = TRUE)))
save(artms_data_randomDF, file = 'data/artms_data_randomDF.RData', 
     compress = 'xz')


## CREATE THE OFFICIAL PHGLOBAL COMING WITH THE PACKAGE-----
setwd("~/sourcecode/artms/ph_full/")

evidence_file <- 'evidence.txt'
keys_file <- 'keys.txt'
contrast_file <- 'contrast.txt'

edf <- read.delim(evidence_file,
                  stringsAsFactors = FALSE)

kdf <- read.delim(keys_file,
                  stringsAsFactors = FALSE,
                  check.names = FALSE)


# Select 2 biological replicates
selectedBR <- c("qx006145", "qx006148", "qx006151", "qx006152")
edfnew <- edf[which(edf$Raw.file %in% selectedBR), ]

# Select columns
edfnew <- edfnew[c("Sequence",
                   "Length",
                   "Modifications",
                   "Modified.sequence",
                   "Oxidation..M..Probabilities",
                   "Phospho..STY..Probabilities",
                   "Oxidation..M..Score.Diffs",
                   "Phospho..STY..Score.Diffs",
                   "Oxidation..M.",
                   "Phospho..STY.",
                   "Missed.cleavages",
                   "Proteins",
                   "Leading.proteins",
                   "Leading.razor.protein",
                   "Type",
                   "Raw.file",
                   "MS.MS.m.z",
                   "Charge",
                   "m.z",
                   "Mass",
                   "Resolution",
                   "Mass.error..ppm.",
                   "Mass.error..Da.",
                   "Uncalibrated.mass.error..ppm.",
                   "Uncalibrated.mass.error..Da.",
                   "Calibrated.retention.time",
                   "Retention.time",
                   "Retention.length",
                   "PEP",
                   "MS.MS.count",
                   "MS.MS.scan.number",
                   "Score",
                   "Delta.score",
                   "Intensity",
                   "Reverse",
                   "Potential.contaminant")]

artms_data_ph_keys <- kdf[which(kdf$RawFile %in% selectedBR),]

# And random sampling lines
n <- round(dim(edfnew)[1] / 25)
artms_data_ph_evidence <- edfnew[sample(nrow(edfnew), n),]

# # print out evidence & keys
# write.table(
#   artms_data_ph_evidence,
#   file = "~/sourcecode/artms/ph/artms_data_ph_evidence.txt",
#   quote = FALSE,
#   sep = "\t",
#   row.names = FALSE,
#   col.names = TRUE
# )
# write.table(
#   artms_data_ph_keys,
#   file = "~/sourcecode/artms/ph/artms_data_ph_keys.txt",
#   quote = FALSE,
#   sep = "\t",
#   row.names = FALSE,
#   col.names = TRUE
# )

# artms_data_ph_evidence----
# artms_data_ph_evidence <- read.delim("~/sourcecode/artms/ph/artms_data_ph_evidence.txt", stringsAsFactors = FALSE)
save(artms_data_ph_evidence, 
     file = '~/github/biodavidjm/artMS/data/artms_data_ph_evidence.RData', 
     compress = TRUE)

# artms_data_ph_keys----
# artms_data_ph_keys <- read.delim("~/sourcecode/artms/extdata/artms_data_ph_keys.txt",
#                                  stringsAsFactors = FALSE)
  
save(artms_data_ph_keys, 
     file = '~/github/biodavidjm/artMS/data/artms_data_ph_keys.RData', 
     compress = TRUE )

# artms_data_ph_contrast -----
contrast_file <- "~/sourcecode/artms/ph_full/contrast.txt"
artms_data_ph_contrast <- readLines(contrast_file, warn = FALSE)
save(artms_data_ph_contrast, 
     file = "~/github/biodavidjm/artMS/data/artms_data_ph_contrast.RData")

# artms_data_ph_config ----
artms_data_ph_config <- artms_config

artms_data_ph_config$files$evidence <- ""
artms_data_ph_config$files$keys <- ""
artms_data_ph_config$files$contrasts <- ""
artms_data_ph_config$files$summary <- ""
artms_data_ph_config$files$output <- "quant-test/results.txt"
artms_data_ph_config$qc$basic <- 0
artms_data_ph_config$qc$extended <- 0
artms_data_ph_config$qc$extendedSummary <- 0
artms_data_ph_config$data$filters$modifications <- "PH"

save(artms_data_ph_config, 
     file = "~/github/biodavidjm/artMS/data/artms_data_ph_config.RData", 
     compress = TRUE)

artms_data_ph_config$files$evidence <- artms_data_ph_evidence
artms_data_ph_config$files$keys <- artms_data_ph_keys
artms_data_ph_config$files$contrasts <- artms_data_ph_contrast
artms_data_ph_config$output_extras <- 0
artms_data_ph_config$msstats$profilePlots <- "before, after"

msresults <- artmsQuantification(yaml_config_file = artms_data_ph_config, 
                                 data_object = TRUE,  
                                 display_msstats = FALSE, 
                                 verbose = TRUE, 
                                 printPDF = FALSE, 
                                 printTables = FALSE)
                            
                    

# Results -----

# Run MSstats on the full version (4 biological replicates)-----
setwd("~/sourcecode/artms/ph_full/")
artmsWriteConfigYamlFile(config_file_name = "artms_full_phglobal.yaml")
artmsQuantification(yaml_config_file = "artms_full_phglobal.yaml")

# Load and write out data objects for artMS
artms_data_ph_msstats_modelqc <- read.delim("~/sourcecode/artms/ph_full/phglobal/new2021/results_ModelQC.txt", stringsAsFactors = FALSE)
save(artms_data_ph_msstats_modelqc, 
     file = "~/github/biodavidjm/artMS/data/artms_data_ph_msstats_modelqc.RData",
     compress = TRUE)

artms_data_ph_msstats_results <-read.delim("~/sourcecode/artms/ph_full/phglobal/new2021/results.txt", stringsAsFactors = FALSE)
save(artms_data_ph_msstats_results, 
     file = "~/github/biodavidjm/artMS/data/artms_data_ph_msstats_results.RData",
     compress = TRUE)


# Using the full version to generate the results file

# To generate results first run:
artms_data_ph_config$files$evidence <- artms_data_ph_evidence
artms_data_ph_config$files$keys <- artms_data_ph_keys
artms_data_ph_config$files$contrasts <- artms_data_ph_contrast

artms_data_ph_quantifications <- artmsQuantification(yaml_config_file = artms_data_ph_config, 
                                                     data_object = TRUE,
                                                     printPDF = TRUE,
                                                     display_msstats = TRUE,
                                                     verbose = TRUE)

artms_data_ph_msstats_results <- as.data.frame(artms_data_ph_quantifications$ComparisonResult)
artms_data_ph_msstats_results$Protein <- as.character(artms_data_ph_msstats_results$Protein)
artms_data_ph_msstats_results$Label <- as.character(artms_data_ph_msstats_results$Label)
artms_data_ph_msstats_results$issue <- as.character(artms_data_ph_msstats_results$issue)

artms_data_ph_msstats_modelqc <- as.data.frame(artms_data_ph_quantifications$ModelQC)
artms_data_ph_msstats_modelqc$RUN <- as.integer(artms_data_ph_msstats_modelqc$RUN)
artms_data_ph_msstats_modelqc$Protein <- as.character(artms_data_ph_msstats_modelqc$Protein)
artms_data_ph_msstats_modelqc$originalRUN <- as.integer(artms_data_ph_msstats_modelqc$originalRUN)
artms_data_ph_msstats_modelqc$GROUP <- as.character(artms_data_ph_msstats_modelqc$GROUP)
artms_data_ph_msstats_modelqc$SUBJECT <- as.character(artms_data_ph_msstats_modelqc$SUBJECT)

save(artms_data_ph_msstats_results, 
     file = "~/github/biodavidjm/artMS/data/artms_data_ph_msstats_results2.RData",
     compress = TRUE)

save(artms_data_ph_msstats_modelqc, 
     file = "~/github/biodavidjm/artMS/data/artms_data_ph_msstats_modelqc2.RData",
     compress = TRUE)



# CORUM dataset----
artms_data_corum_mito_database <- read.delim("~/github/biodavidjm/artMS/inst/extdata/20170801_corum_mitoT.txt", 
                                             stringsAsFactors = FALSE)
  
save(artms_data_corum_mito_database,
     file = 'data/artms_data_corum_mito_database.RData',
     compress = TRUE)





# PATHOGENS----
message("--- PATHOGEN IN SAMPLES: TB\n")
artms_data_pathogen_TB <- read.delim('~/Box Sync/db/uniprot/uniprot-tr-myctb_tuberculosis_ATCC35801_TMC10-onlyEntryID.fasta',
                                     header = FALSE,
                                     sep = "\t",
                                     quote = "",
                                     stringsAsFactors = FALSE) # pathogen.ids$Entry, "TB",
  
names(artms_data_pathogen_TB) <- c('Entry')
save(artms_data_pathogen_TB, 
     file = '~/github/biodavidjm/artMS/data/artms_data_pathogen_TB.RData', 
     compress = 'xz')

message("--- PATHOGEN IN SAMPLES: LEGIONELLA PNEUMOPHILA")
artms_data_pathogen_LPN <-
  read.delim(
    '~/Box Sync/db/uniprot/uniprot-legionella-proteome_UP000000609.txt',
    header = TRUE,
    sep = "\t",
    quote = "",
    stringsAsFactors = FALSE
  ) # pathogen.ids$Entry, "Lpn",
artms_data_pathogen_LPN <- artms_data_pathogen_LPN[c('Entry')]

save(artms_data_pathogen_LPN,
     file = '~/github/biodavidjm/artMS/data/artms_data_pathogen_LPN.RData',
     compress = 'xz')


# VIGNETTES
artmsQualityControlEvidenceBasic(evidence_file = artms_data_ph_evidence,
                                 keys_file = artms_data_ph_keys,
                                 prot_exp = "PH")







