# CREATE SAMPLE DATA FILES FOR THE PACKAGE

library(biomaRt)
library(bit64)
library(circlize)
library(cluster)
library(corrplot)
library(data.table)
library(factoextra)
library(FactoMineR)
library(getopt)
library(ggdendro)
library(ggplot2)
library(gProfileR)
library(grid)
library(limma)
library(MSstats)
library(openxlsx)
library(org.Ag.eg.db)
library(org.At.tair.db)
library(org.Bt.eg.db)
library(org.Ce.eg.db)
library(org.Cf.eg.db)
library(org.Dm.eg.db)
library(org.Dr.eg.db)
library(org.EcK12.eg.db)
library(org.EcSakai.eg.db)
library(org.Gg.eg.db)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Mmu.eg.db)
library(org.Pf.plasmo.db)
library(org.Pt.eg.db)
library(org.Rn.eg.db)
library(org.Ss.eg.db)
library(org.Xl.eg.db)
library(PerformanceAnalytics)
library(pheatmap)
library(plotly)
library(plyr)
library(RColorBrewer)
library(reshape2)
library(seqinr)
library(shiny)
library(stringr)
library(VennDiagram)
library(yaml)
library(graphics)
library(grDevices)
library(stats)
library(utils)

setwd('~/github/biodavidjm/artMS/')

# GENERATE RANDOM FILE
randomDF <- data.frame(replicate(10,sample(0:1,100,rep=TRUE)))
save(randomDF, file = 'data/randomDF.RData', compress = 'xz')

# PH FILES
# MaxQuant Evidence file
ph_evidence <- read.delim('inst/extdata/evidence.txt', stringsAsFactors = F)
save(ph_keys, file='data/ph_keys.RData', compress = 'xz')

# Keys file (experimental design)
ph_keys <- read.delim("inst/extdata/keys.txt", stringsAsFactors = F)
save(ph_evidence, file = 'data/ph_evidence.RData', compress = 'xz')

# CORUM dataset
corum_mito_database <- read.delim("inst/extdata/20170801_corum_mitoT.txt", stringsAsFactors = F)
save(corum_mito_database, file = 'data/corum_mito_database.RData', compress = 'xz')

# CONFIGURATION FILE
artms_config <- yaml.load_file("inst/extdata/artms_config.yaml")
save(artms_config, file = 'data/artms_config.RData', compress = 'xz')

load("data/artms_config.RData")



# Testing artMS

## FRACTIONS
evidence_file <- '~/experiments/artms/fractions/petroski-cul4-evidence.txt'
keys_file <- '~/experiments/artms/fractions/petroski-cul4-keys.txt'
contrast_file <- '~/experiments/artms/fractions/petroski-cul4-contrast.txt'
yaml_config_file <- '~/experiments/artms/fractions/results/ab20180402/config-petroski-debugging.yaml'

## PH
evidence_file <- '~/experiments/artms/ph/evidence.txt'
keys_file <- '~/experiments/artms/ph/keys.txt'
contrast_file <- '~/experiments/artms/ph/contrast.txt'
## PHSITES
yaml_config_file <- '~/experiments/artms/ph/phsites/phsites_config.yaml'
artms_main(yaml_config_file = yaml_config_file)

setwd('~/experiments/artms/ph/phsites/')
log2fc_file = "phsites-results.txt"
modelqc_file = "phsites-results_ModelQC.txt"
specie = "human"
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
evidence_file <- '~/experiments/artms/silac/RI__Endosome_Abundance_NoDrugvsDrug-evidence.txt'
keys_file <- '~/experiments/artms/silac/RI__Endosome_Abundance_NoDrugvsDrug-keys.txt'
contrast_file <- '~/experiments/artms/silac/RI__Endosome_Abundance_NoDrugvsDrug-contrasts.txt'
yaml_config_file <- '~/experiments/artms/silac/results/config-silac.yaml'

## ABUNDANCE, technical replicates
evidence_file <- '~/experiments/artms/technical_replicas/201706-FLU-HTBE-H5N1-AB-evidence.txt'
keys_file <- '~/experiments/artms/technical_replicas/FLU-HTBE-H5N1-AB-keys.txt'
contrast_file <- '~/experiments/artms/technical_replicas/FLU-HTBE-H5N1-contrasts-final.txt'
yaml_config_file <- '~/experiments/artms/technical_replicas/configTR.yaml'

# Abundance, no technical replicates
evidence_file <- '~/experiments/artms/thp1_ab_h1n1/FLU-THP1-H1N1-AB-evidence.txt'
keys_file <- '~/experiments/artms/thp1_ab_h1n1/FLU-THP1-H1N1-AB-keys.txt'
contrast_file <- '~/experiments/artms/thp1_ab_h1n1/FLU-THP1-H1N1-AB-contrasts.txt'

## ANALYSIS OF QUANTIFICATIONS
setwd('~/experiments/artms/thp1_ab_h1n1/results/testing/')

artms_analysisQuantifications(log2fc_file = "ab-testing-new-results.txt",
                              modelqc_file = "ab-testing-new-results_ModelQC.txt",
                              specie = "human",
                              isPtm = "noptm",
                              enrich = TRUE,
                              output_dir = "AnalysisQuantifications",
                              isFluomics = TRUE,
                              isBackground = "nobackground",
                              mnbr = 2,
                              l2fc_thres = 1,
                              ipval = "pvalue",
                              pathogen = "nopathogen")

log2fc_file = "ab-testing-new-results.txt"
modelqc_file = "ab-testing-new-results_ModelQC.txt"
specie = "human"
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
  input_file = '~/experiments/artms/technical_replicas/results/FLU-HTBE-H5N1-results.txt', 
  output_file = '~/experiments/artms/technical_replicas/results/FLU-HTBE-H5N1-results-plotheatmap.pdf')

## Evidence to MIST and MISTIN
artms_evidenceToMISTformat(metric = "int", 
                           input_file = '~/experiments/artms/technical_replicas/201706-FLU-HTBE-H5N1-AB-evidence.txt', 
                           output_file = '~/experiments/artms/technical_replicas/201706-FLU-HTBE-H5N1-AB-evidence-mist-int.txt', 
                           keys_file = '~/experiments/artms/technical_replicas/FLU-HTBE-H5N1-AB-keys.txt', 
                           species = 'HUMAN-FLUOMICS', 
                           uniprot_dir = '~/Box Sync/db/mist/')

## Evidence QC
artms_evidenceQC(evidence_file = evidence_file, 
                 keys_file = keys_file, 
                 prot_exp = "ph")




## ANNOTATIONS
# Generate the annotation system
symbols <- c('JAK1','AATK','A2BP1','A2LD1')

select(org.Hs.eg.db, symbols, c("ENTREZID","GENENAME"), "ALIAS")

# RANDOMLY SELECT KEYS FROM UNIPROT HUMANS
uniprots <- c("Q6P996")
uniprots <- Rkeys(org.Hs.egUNIPROT)[1:100]

ano <- artms_mapUniprot2entrezGeneName(theUniprots = uniprots, specie = "human")

library(org.Hs.eg.db)
library(org.Mm.eg.db)

# UNIPROT TO ENTREZ
uni2entrez <- select(org.Hs.eg.db, uniprots, "ENTREZID", "UNIPROT")

# Uniprot 2 symbol
mappings <- select(org.Hs.eg.db, uniprots, c("UNIPROT", "SYMBOL", "GENENAME"), keytype = "UNIPROT")
mappings <- AnnotationDbi::select(org.Hs.eg.db, uniprots, c("UNIPROT", "SYMBOL", "GENENAME"), keytype = "UNIPROT")
# Remove redundancies
mappings <- mappings[!duplicated(mappings$UNIPROT),]


# IF I NEED TO BRING ALL THE UNIPROT IDS
uni_ids <- keys(org.Hs.eg.db, c("UNIPROT"))



