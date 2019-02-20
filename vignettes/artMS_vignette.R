## ---- echo = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  eval=FALSE
)

## ----library---------------------------------------------------------------
#  library(artMS)

## ---- eval = TRUE----------------------------------------------------------
library(artMS)
artmsWriteConfigYamlFile(config_file_name = "config.yaml", 
                         verbose = FALSE)

## ---- eval = TRUE----------------------------------------------------------
# But for illustration purposes printing only INTDIST plot:
library(artMS)
suppressWarnings(
artmsQualityControlEvidenceBasic(evidence_file = artms_data_ph_evidence,
                                 keys_file = artms_data_ph_keys,
                                 prot_exp = "PH",
                                 plotINTDIST = TRUE,
                                 plotREPRO = FALSE,
                                 plotCORMAT = FALSE,
                                 plotINTMISC = FALSE,
                                 plotPTMSTATS = FALSE,
                                 printPDF = FALSE,
                                 verbose = FALSE))


## ---- eval = FALSE---------------------------------------------------------
#  # This example adds annotations to the evidence file available in
#  # artMS, based on the column 'Proteins'.
#  
#  evidence_anno <- artmsAnnotationUniprot(x = artms_data_ph_evidence,
#                                          columnid = 'Proteins',
#                                          species = 'human')

## ---- eval = FALSE---------------------------------------------------------
#  artms_data_ph_evidence <- artmsChangeColumnName(
#                                 dataset = artms_data_ph_evidence,
#                                 oldname = "Phospho..STY.",
#                                 newname = "PH_STY")

## ---- eval=FALSE-----------------------------------------------------------
#  # The data must be annotated (Protein and Gene columns)
#  data_annotated <- artmsAnnotationUniprot(
#                        x = artms_data_ph_msstats_results,
#                        columnid = "Protein",
#                        species = "human")
#  # And then the enrichment
#  enrich_set <- artmsEnrichLog2fc(
#                     dataset = data_annotated,
#                     species = "human",
#                     background = unique(data_annotated$Gene),
#                     verbose = FALSE)

## ---- eval=FALSE-----------------------------------------------------------
#  # annotate the MSstats results to get the Gene name
#  data_annotated <- artmsAnnotationUniprot(
#                                       x = artms_data_ph_msstats_results,
#                                       columnid = "Protein",
#                                       species = "human")
#  
#  # Filter the list of genes with a log2fc > 2
#  filtered_data <-
#  unique(data_annotated$Gene[which(data_annotated$log2FC > 2)])
#  
#  # And perform enrichment analysis
#  data_annotated_enrich <- artmsEnrichProfiler(
#                                     x = filtered_data,
#                                     categorySource = c('KEGG'),
#                                     species = "hsapiens",
#                                     background = unique(data_annotated$Gene))

