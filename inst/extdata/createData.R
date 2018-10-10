# CREATE SAMPLE DATA FILES FOR THE PACKAGE

# library(biomaRt)
# library(bit64)
# library(circlize)
# library(cluster)
# library(corrplot)
# library(data.table)
# library(factoextra)
# library(FactoMineR)
# library(getopt)
# library(ggdendro)
# library(ggplot2)
# library(gplots)
# library(ggrepel)
# library(gProfileR)
# library(grid)
# library(limma)
# library(MSstats)
# library(openxlsx)
# library(org.Ag.eg.db)
# library(org.At.tair.db)
# library(org.Bt.eg.db)
# library(org.Ce.eg.db)
# library(org.Cf.eg.db)
# library(org.Dm.eg.db)
# library(org.Dr.eg.db)
# library(org.EcK12.eg.db)
# library(org.EcSakai.eg.db)
# library(org.Gg.eg.db)
# library(org.Hs.eg.db)
# library(org.Mm.eg.db)
# library(org.Mmu.eg.db)
# library(org.Pf.plasmo.db)
# library(org.Pt.eg.db)
# library(org.Rn.eg.db)
# library(org.Sc.sgd.db)
# library(org.Ss.eg.db)
# library(org.Xl.eg.db)
# library(PerformanceAnalytics)
# library(pheatmap)
# library(plotly)
# library(plyr)
# library(RColorBrewer)
# library(reshape2)
# library(seqinr)
# library(stringr)
# library(VennDiagram)
# library(yaml)
# library(graphics)
# library(grDevices)
# library(stats)
# library(UpSetR)
# library(utils)
# library(formatR)


# TESTING OPTIONS


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE DATA

setwd('~/github/biodavidjm/artMS/')

# GENERATE RANDOM FILE
artms_data_randomDF <- data.frame(replicate(10, sample(0:1, 100, rep = TRUE)))
save(artms_data_randomDF, file = 'data/artms_data_randomDF.RData', 
     compress = 'xz')

# PH FILES

# Reduced version of an Evidence file (generated below)
artms_data_ph_evidence <- read.delim(
  "~/sourcecode/artms/ph/artms_data_ph_evidence.txt",
  stringsAsFactors = FALSE)
save(artms_data_ph_evidence, file = 'data/artms_data_ph_evidence.RData', 
     compress = 'xz')

# Reduced version of the Keys file (experimental design)
artms_data_ph_keys <-
  read.delim("~/sourcecode/artms/extdata/artms_data_ph_keys.txt",
             stringsAsFactors = FALSE)
save(artms_data_ph_keys, file = 'data/artms_data_ph_keys.RData', 
     compress = 'xz')

# Reduced version of the results
artms_data_ph_msstats_results <-
  read.delim(
    "~/sourcecode/artms/extdata/artms_data_ph_msstats_results.txt",
    stringsAsFactors = FALSE
  )
save(artms_data_ph_msstats_results,
     file = 'data/artms_data_ph_msstats_results.RData',
     compress = 'xz')

artms_data_ph_proteinGroups <-
  read.delim("~/sourcecode/artms/ph/proteinGroups.txt",
             stringsAsFactors = FALSE)
save(artms_data_ph_proteinGroups,
     file = 'data/artms_data_ph_proteinGroups.RData',
     compress = 'xz')

# CORUM dataset
artms_data_corum_mito_database <-
  read.delim("inst/extdata/20170801_corum_mitoT.txt", stringsAsFactors = FALSE)
save(artms_data_corum_mito_database,
     file = 'data/artms_data_corum_mito_database.RData',
     compress = 'xz')

# CONFIGURATION FILE
# library(yaml)
artms_config <- yaml.load_file("inst/extdata/artms_config.yaml")
save(artms_config, file = 'data/artms_config.RData', compress = 'xz')

# PATHOGENS

cat("--- PATHOGEN IN SAMPLES: TB\n")
artms_data_pathogen_TB <-
  read.delim(
    '~/Box Sync/db/uniprot/uniprot-tr-myctb_tuberculosis_ATCC35801_TMC10-onlyEntryID.fasta',
    header = FALSE,
    sep = "\t",
    quote = "",
    stringsAsFactors = FALSE
  ) # pathogen.ids$Entry, "TB",
names(artms_data_pathogen_TB) <- c('Entry')
save(artms_data_pathogen_TB, 
     file = '~/github/biodavidjm/artMS/data/artms_data_pathogen_TB.RData', 
     compress = 'xz')

cat("--- PATHOGEN IN SAMPLES: LEGIONELLA PNEUMOPHILA\n")
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Warning testing zone

x <- list(a = 1, b = 1:3, c = 10:100)
vapply(x, FUN = length, FUN.VALUE = 0)

sapply(x, FUN = length)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# VIGNETTES

artms_qualityControlEvidenceBasic(evidence_file = artms_data_ph_evidence,
                                  keys_file = artms_data_ph_keys,
                                  prot_exp = "PH")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TESTING artMS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# QC PLOTS EXTENDED
setwd("~/Box Sync/tempStuff/david2alex/alex2david/1D/")
evidence_file <- "evidence.txt"
keys_file <- "keys.txt"
summary_file <- "summary.txt"

artms_qualityControlEvidenceExtended(evidence_file = "evidence.txt",
                                     keys_file = "keys.txt")

artms_qualityControlSummaryExtended(summary_file = "summary.txt",
                                    keys_file = "keys.txt")

artms_quantification(yaml_config_file = "config.yaml")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## APMS FLUOMICS
setwd("~/sourcecode/artms/apms/")

here <- read.delim("a549-PB1-keys.txt", stringsAsFactors = FALSE)

artms_qualityControlEvidenceBasic(evidence_file = "a549-PB1-evidence.txt",
                                  keys_file = "a549-PB1-keys.txt",
                                  prot_exp = "APMS")

artms_qualityControlEvidenceExtended(evidence_file = "a549-PB1-evidence.txt",
                                     keys_file = "a549-PB1-keys.txt")

artms_qualityControlSummaryExtended(summary_file = "summary.txt",
                                    keys_file = "a549-PB1-keys.txt")

artms_quantification(yaml_config_file = "apms_config.yaml")

artms_evidenceToSAINTqFormat(evidence_file = "a549-PB1-evidence.txt",
                             keys_file = "a549-PB1-keys.txt",
                             output_dir = "saintq_folder")

artms_evidenceToSaintExpressFormat(
  evidence_file = "a549-PB1-evidence.txt",
  keys_file = "a549-PB1-keys.txt",
  output_file = "a549-PB1-saintexpress.txt",
  ref_proteome_file = "~/Box Sync/db/flu/fluomics-uniprot-hsa_20170516.fasta"
)

artms_msstats_summary(
  evidence_file = "a549-PB1-evidence.txt",
  keys_file = "a549-PB1-keys.txt",
  prot_group_file = "proteinGroups.txt"
)

setwd("~/sourcecode/artms/apms/results/")
artms_analysisQuantifications(
  # log2fc_file = "a549-PB1-results.txt",
  # modelqc_file = "a549-PB1-results_ModelQC.txt",
  species = "HUMAN",
  output_dir = "analysis"
)

artms_volcanoPlot(
  mss_results = "a549-PB1-results.txt",
  lfc_upper = 1,
  lfc_lower = -1,
  FDR = 0.05,
  output_name = "a549-PB1-results-volcanoPlot.pdf",
  PDF = TRUE
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## FRACTIONS
setwd('~/sourcecode/artms/fractions/')
evidence_file <-
  '~/sourcecode/artms/fractions/petroski-cul4-evidence.txt'
keys_file <- '~/sourcecode/artms/fractions/petroski-cul4-keys.txt'
contrast_file <-
  '~/sourcecode/artms/fractions/petroski-cul4-contrast.txt'

yaml_config_file <-
  '~/sourcecode/artms/fractions/config-petroski-debugging.yaml'

# Quantifications
artms_quantification(yaml_config_file = yaml_config_file)

# Analysis of Quantifications

setwd('~/sourcecode/artms/fractions/results/')
artms_analysisQuantifications(
  log2fc_file = "petroski-cul4-debug2-results.txt",
  modelqc_file = "petroski-cul4-debug2-results_ModelQC.txt",
  species = "human",
  output_dir = "analysis"
)

setwd('~/sourcecode/artms/fractions/')
artms_msstats_summary(
  evidence_file = "petroski-cul4-evidence.txt",
  keys_file = "petroski-cul4-keys.txt",
  prot_group_file = "proteinGroups.txt",
  results_file = "results/petroski-cul4-debug2-results.txt"
)


#-------------------------------------------------------------------------------
## CREATE THE OFFICIAL PHGLOBAL COMING WITH THE PACKAGE
evidence_file <- 'evidence.txt'
keys_file <- 'keys.txt'
contrast_file <- 'contrast.txt'

edf <-
  read.delim(evidence_file,
             stringsAsFactors = FALSE,
             check.names = FALSE)
kdf <-
  read.delim(keys_file,
             stringsAsFactors = FALSE,
             check.names = FALSE)

# Select 2 biological replicates
selectedBR <- c("qx006145", "qx006148", "qx006151", "qx006152")
edfnew <- edf[which(edf$`Raw file` %in% selectedBR), ]
kdfnew <- kdf[which(kdf$RawFile %in% selectedBR),]

# And random sampling lines
n <- round(dim(edfnew)[1] / 7)
edfnew2 <- edfnew[sample(nrow(edfnew), n),]

# print out evidence & keys
write.table(
  edfnew2,
  file = "~/sourcecode/artms/ph/artms_data_ph_evidence.txt",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)
write.table(
  kdfnew,
  file = "~/github/biodavidjm/artMS/inst/extdata/artms_data_ph_keys.txt",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)

#-------------------------------------------------------------------------------
# PH GLOBAL:
setwd('~/sourcecode/artms/ph/')
contrast_file <- 'contrast.txt'
evidence_file = "evidence.txt"
prot_group_file = "proteinGroups.txt"
keys_file = "keys.txt"
contrast_file <- 'contrast.txt'

artms_qualityControlSummaryExtended(summary_file = "summary.txt", 
                                    keys_file = "keys.txt")

artms_writeConfigYamlFile()

artms_quantification("phglobal/phglobal_config.yaml")

artms_msstats_summary(
  evidence_file = "evidence.txt",
  prot_group_file = "proteinGroups.txt",
  keys_file = "keys.txt",
  results_file = "phglobal/phglobal-results.txt"
)

artms_replicatePlots(
  input_file = evidence_file,
  keys_file = keys_file,
  replicate_file = "replicates_plots.txt",
  prot_exp = "PH",
  out_file = "replicate-plots.txt"
)

setwd('~/sourcecode/artms/ph/phglobal/')

artms_analysisQuantifications(
  log2fc_file = "phglobal-results.txt",
  modelqc_file = "phglobal-results_ModelQC.txt",
  species = "human",
  output_dir = "analysis_name"
)

#-------------------------------------------------------------------------------
# PH SITES
setwd('~/sourcecode/artms/ph/')

artms_proteinToSiteConversion(
  evidence_file = "evidence.txt",
  ref_proteome_file = "uniprot_canonical.fasta",
  output_file = "evidence-sites.txt",
  mod_type = "ph"
)

artms_qualityControlEvidenceExtended(evidence_file = "evidence-sites.txt",
                                     keys_file = "keys.txt")

artms_quantification(yaml_config_file = "phsites/phsites_config.yaml")

setwd('~/sourcecode/artms/ph/phsites/')
artms_analysisQuantifications(
  log2fc_file = "phsites-results.txt",
  modelqc_file = "phsites-results_ModelQC.txt",
  species = "human",
  output_dir = "analysisQuant2",
  isPtm = "ptmsites",
  enrich = FALSE
)

log2fc_file = "phsites-results.txt"
modelqc_file = "phsites-results_ModelQC.txt"
species = "human"
output_dir = "analysisQuant"
isPtm = "ptmsites"
enrich = FALSE
l2fc_thres = 1.5
choosePvalue = "adjpvalue"
isBackground = "nobackground"
mnbr = 2
isFluomics = FALSE
pathogen = "nopathogen"

df = imputedDF
pathogen = pathogen
species = species
ptmType = isPtm
output_name = log2fc_file


setwd('~/sourcecode/artms/ph/phsites/analysisQuant2_adjpvalue/')

filename <- "phsites-results-imputedL2fcExtended.txt"

artmsPhosfateOutput(inputFile = filename)
artmsPhotonOutput(inputFile = filename)

here <- artms_generatePhSiteExtended(df = "phsites-results-abundance-long.txt", 
                             pathogen = "nopathogen", 
                             species = "human", 
                             ptmType = "ptmsites", 
                             output_name = "whatever")

setwd("/Users/djm75/Box Sync/projects/FluomicsProteomics/Flu-human-exvivo/PTMs/HTBE/H1N1/ph_2017/results/20171015sites/a20180205_pvalue")
filename <- "201710-FLU-HTBE-H1N1-PH-modacc-results-imputedL2fcExtended.txt"
artmsPhosfateOutput(inputFile = filename)

#-------------------------------------------------------------------------------
# PH REDUCED
setwd('~/sourcecode/artms/ph/phglobalreduced/')

artms_quantification("phglobal_reduced_config.yaml")

artms_analysisQuantifications(
  log2fc_file = "ph-reduced-results.txt",
  modelqc_file = "ph-reduced-results_ModelQC.txt",
  species = "human",
  isPtm = "global",
  enrich = TRUE,
  output_dir = "testingARTMS3",
  mnbr = 2,
  l2fc_thres = 1.5,
  ipval = "pvalue"
)

#-------------------------------------------------------------------------------
ph_results_wide <- artms_resultsWide(
  results_msstats = artms_data_ph_msstats_results,
  output_file = NULL)

summary_spectral_counts <-
  artms_spectralCounts(evidence_file = artms_data_ph_evidence,
                       keys_file = artms_data_ph_keys)

evidence_anno <-
  artms_annotationUniprot(data = artms_data_ph_evidence,
                          columnid = "Proteins",
                          sps = "human")

uniprots_anno <- artms_mapUniprot2entrezGeneName(
  uniprotkb = unique(artms_data_ph_evidence$Proteins),
  species = "human")

data_annotated <-
  artms_annotationUniprot(data = artms_data_ph_msstats_results,
                          columnid = "Protein",
                          sps = "human")

# Filter the list of genes with a log2fc > 2
filtered_data <-
  unique(data_annotated$Gene[which(data_annotated$log2FC > 2)])

# And enrich it
data_annotated_enrich <- artms_enrichProfiler(
  x = filtered_data,
  categorySource = c('KEGG'),
  species = "hsapiens",
  background = unique(data_annotated$Gene)
)

# -----------------------------------------------------------------------------
artms_plotHeatmapQuant(
  input_file = artms_data_ph_msstats_results,
  species = "human",
  output_file = NULL,
  whatPvalue = "pvalue",
  lfc_lower = -1,
  lfc_upper = 1
)
# -----------------------------------------------------------------------------
artms_volcanoPlot(mss_results = artms_data_ph_msstats_results,
                  whatPvalue = "pvalue",
                  PDF = FALSE)

# -----------------------------------------------------------------------------
# The data must be annotated (Protein and Gene columns)
data_annotated <- artms_annotationUniprot(data = artms_data_ph_msstats_results,
                                          columnid = "Protein",
                                          sps = "human")
# And then the enrichment
enrich_set <- artms_enrichLog2fc(
  dataset = data_annotated,
  species = "human",
  background = unique(data_annotated$Gene),
  heatmaps = TRUE
)

dataset = data_annotated
species = "human"
background = unique(data_annotated$Gene)
heatmaps = TRUE

#----------------------------------------------------------------------
artms_isEvidenceNewVersion(evidence_file = artms_data_ph_evidence)


#----------------------------------------------------------------------
# Adding a new column with the main species of the data. Easy.
# But the main functionality is to add both the host-species and a pathogen,
# which is not illustrated in this example
artms_annotateSpecie(df = artms_data_ph_msstats_results, species = "human")

#-------------------------------------------------------------------------------
# First, let's make the "replicate file" (in a data.frame)
x_names <-
  c("condition1",
    "rep1_1",
    "rep1_2",
    "condition2",
    "rep2_1",
    "rep2_2")
x_values <-
  c("Cal33", "Cal33-1", "Cal33-4", "HSC6", "HSC6-2", "HSC6-3")
replica_info <- data.frame(t(x_values))
colnames(replica_info) <- x_names

# Now let's make the plots (it is recommended to use the <out_file> option
# and print the results to a file)
artms_replicatePlots(
  input_file = artms_data_ph_evidence,
  keys_file = artms_data_ph_keys,
  replicate_file = replica_info,
  out_file = "whatever.txt",
  prot_exp = "PH"
)



## PHSITES
setwd('~/sourcecode/artms/ph/')
artms_proteinToSiteConversion(
  evidence_file = "evidence.txt",
  ref_proteome_file = "uniprot_canonical.fasta",
  output_file = "phsite_evidence.txt",
  mod_type = "ph"
)


# Generate the site-evidence.txt file:
artms_proteinToSiteConversion(evidence_file = '')

yaml_config_file <-
  '~/sourcecode/artms/ph/phsites/phsites_config.yaml'
artms_main(yaml_config_file = yaml_config_file)

setwd('~/sourcecode/artms/ph/phsites/')
log2fc_file = "phsites-results.txt"
modelqc_file = "phsites-results_ModelQC.txt"
species = "human"
enrich = "yesenrich"
output_dir = "resultsTesting"
isFluomics = TRUE
isPtm = "yesptmph"
isBackground = "nobackground"
mnbr = 2
threshold = 1
ipval = "pvalue"
pathogen = "nopathogen"

## SILAC
evidence_file <-
  '~/sourcecode/artms/silac/RI__Endosome_Abundance_NoDrugvsDrug-evidence.txt'
keys_file <-
  '~/sourcecode/artms/silac/RI__Endosome_Abundance_NoDrugvsDrug-keys.txt'
contrast_file <-
  '~/sourcecode/artms/silac/RI__Endosome_Abundance_NoDrugvsDrug-contrasts.txt'
yaml_config_file <-
  '~/sourcecode/artms/silac/results/config-silac.yaml'

evidence2silac <-
  artms_SILACtoLong(evidence_file = evidence_file, 
                    output = "silac-evidence.txt")

## ABUNDANCE, technical replicates
evidence_file <-
  '~/sourcecode/artms/technical_replicas/201706-FLU-HTBE-H5N1-AB-evidence.txt'
keys_file <-
  '~/sourcecode/artms/technical_replicas/FLU-HTBE-H5N1-AB-keys.txt'
contrast_file <-
  '~/sourcecode/artms/technical_replicas/FLU-HTBE-H5N1-contrasts-final.txt'
yaml_config_file <-
  '~/sourcecode/artms/technical_replicas/configTR.yaml'

# Abundance, no technical replicates
evidence_file <-
  '~/sourcecode/artms/thp1_ab_h1n1/FLU-THP1-H1N1-AB-evidence.txt'
keys_file <-
  '~/sourcecode/artms/thp1_ab_h1n1/FLU-THP1-H1N1-AB-keys.txt'
contrast_file <-
  '~/sourcecode/artms/thp1_ab_h1n1/FLU-THP1-H1N1-AB-contrasts.txt'

## ANALYSIS OF QUANTIFICATIONS
setwd('~/sourcecode/artms/thp1_ab_h1n1/results/testing/')

artms_analysisQuantifications(
  log2fc_file = "ab-testing-new-results.txt",
  modelqc_file = "ab-testing-new-results_ModelQC.txt",
  species = "human",
  isPtm = "noptm",
  enrich = TRUE,
  output_dir = "AnalysisQuantifications",
  isFluomics = TRUE,
  isBackground = "nobackground",
  mnbr = 2,
  l2fc_thres = 1,
  ipval = "pvalue",
  pathogen = "nopathogen"
)

q <- resultsHeatmap(results_file = "ab-testing-new-results.txt",
                    save_file = "whatever.pdf",
                    species = "human")
print(q)

artms_resultsWide(evidence_file = "results/testing/ab-testing-new-results.txt",
              output_file = "results/testing/ab-testing-new-results-wide.txt")

artms_dataPlots(
  input_file = "results/testing/ab-testing-new-results-mss-normalized.txt",
  output_file = "results/testing/ab-testing-new-results-mss-normalized.pdf")

artms_plotHeatmapQuant(input_file = "ab-testing-new-results.txt",
                       species = "human")

print(here)
artms_msstats_summary(
  evidence_file = "FLU-THP1-H1N1-AB-evidence.txt",
  prot_group_file = "proteinGroups.txt",
  keys_file = "FLU-THP1-H1N1-AB-keys.txt",
  results_file = "results/testing/ab-testing-new-results.txt",
  return_df = TRUE
)

artms_spectralCounts(evidence_file = "FLU-THP1-H1N1-AB-evidence.txt",
                     keys_file = "FLU-THP1-H1N1-AB-keys.txt",
                     output_file = "FLU-THP1-H1N1-AB-spectral_counts.txt")

evidence <-
  read.delim("FLU-THP1-H1N1-AB-evidence.txt", stringsAsFactors = FALSE)
keys <-
  read.delim("FLU-THP1-H1N1-AB-keys.txt", stringsAsFactors = FALSE)
evidenceKeys <-
  artms_mergeEvidenceAndKeys(data = evidence, keys = keys)

evidenceKeys <-
  artms_mergeEvidenceKeysByFiles(
    evidence_file = "FLU-THP1-H1N1-AB-evidence.txt", 
    keys_file = "FLU-THP1-H1N1-AB-keys.txt")

evidence_filtered <- artms_filterMaxqData(data = evidence)


evidence_file = "FLU-THP1-H1N1-AB-evidence.txt"
prot_group_file = "proteinGroups.txt"
keys_file = "FLU-THP1-H1N1-AB-keys.txt"
results_file = "results/testing/ab-testing-new-results.txt"
return_results = TRUE


log2fc_file = "ab-testing-new-results.txt"
modelqc_file = "ab-testing-new-results_ModelQC.txt"
species = "human"
isPtm = "noptm"
enrich = TRUE
output_dir = "AnalysisQuantifications"
isFluomics = TRUE
isBackground = "nobackground"
mnbr = 2
l2fc_thres = 1
ipval = "pvalue"
pathogen = "nopathogen"


## CONFIG LOADING
# yaml_config_file <- '~/github/biodavidjm/artMS/data-raw/artms_config.yaml'
artms_main(yaml_config_file)

## Testing individual functions
here <- artms_plotHeatmap(
  input_file = '~/sourcecode/artms/technical_replicas/results/FLU-HTBE-H5N1-results.txt',
  output_file = '~/sourcecode/artms/technical_replicas/results/FLU-HTBE-H5N1-results-plotheatmap.pdf')

## Evidence to MIST and MISTIN
artms_evidenceToMISTformat(
  quant_variable = "int",
  input_file = '~/sourcecode/artms/technical_replicas/201706-FLU-HTBE-H5N1-AB-evidence.txt',
  output_file = '~/sourcecode/artms/technical_replicas/201706-FLU-HTBE-H5N1-AB-evidence-mist-int.txt',
  keys_file = '~/sourcecode/artms/technical_replicas/FLU-HTBE-H5N1-AB-keys.txt',
  species = 'HUMAN-FLUOMICS',
  uniprot_dir = '~/Box Sync/db/mist/'
)

## Evidence QC
artms_qualityControlEvidenceBasic(evidence_file = evidence_file,
                                  keys_file = keys_file,
                                  prot_exp = "ph")




## ANNOTATIONS
# Generate the annotation system
symbols <- c('JAK1', 'AATK', 'A2BP1', 'A2LD1')

select(org.Hs.eg.db, symbols, c("ENTREZID", "GENENAME"), "ALIAS")

# RANDOMLY SELECT KEYS FROM UNIPROT HUMANS
uniprots <- as.list(Rkeys(org.Hs.egUNIPROT)[5000:5050])

exampleID <- c("Q6P996", "B1N8M6")
artmsMapUniprot2Entrez(uniprotkb = exampleID, 
                                          species = "HUMAN")
df_example_anno

ano <-
  artms_mapUniprot2entrezGeneName(uniprotkb = uniprots, species = "human")

# library(org.Hs.eg.db)
# library(org.Mm.eg.db)

# UNIPROT TO ENTREZ
uni2entrez <- select(org.Hs.eg.db, uniprots, "ENTREZID", "UNIPROT")

# Uniprot 2 symbol
mappings <-
  select(org.Hs.eg.db,
         uniprots,
         c("UNIPROT", "SYMBOL", "GENENAME"),
         keytype = "UNIPROT")
mappings <-
  AnnotationDbi::select(org.Hs.eg.db,
                        uniprots,
                        c("UNIPROT", "SYMBOL", "GENENAME", "ENTREZID"),
                        keytype = "UNIPROT")
# Remove redundancies
mappings <- mappings[!duplicated(mappings$UNIPROT), ]


# IF I NEED TO BRING ALL THE UNIPROT IDS
uni_ids <- keys(org.Hs.eg.db, c("UNIPROT"))
