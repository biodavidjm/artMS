## ---- echo = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  eval=FALSE
)

## ---- eval = FALSE---------------------------------------------------------
#  getRversion()

## ---- eval = FALSE---------------------------------------------------------
#  source("https://bioconductor.org/biocLite.R")
#  biocLite(c('org.Ag.eg.db', 'org.At.tair.db', 'org.Bt.eg.db',
#             'org.Ce.eg.db', 'org.Cf.eg.db', 'org.Dm.eg.db',
#             'org.Dr.eg.db', 'org.EcK12.eg.db', 'org.EcSakai.eg.db',
#             'org.Gg.eg.db', 'org.Hs.eg.db', 'org.Mm.eg.db',
#             'org.Mmu.eg.db', 'org.Pf.plasmo.db', 'org.Pt.eg.db',
#             'org.Rn.eg.db', 'org.Sc.sgd.db', 'org.Ss.eg.db',
#             'org.Xl.eg.db'))

## ---- eval = FALSE---------------------------------------------------------
#  install.packages("devtools")
#  library(devtools)
#  install_github("biodavidjm/artMS", build_vignettes=TRUE)

## ---- eval = FALSE---------------------------------------------------------
#  library(artMS)
#  ?artms_qualityControlEvidenceBasic

## ---- eval = FALSE---------------------------------------------------------
#  # First go to a local working directory: several pdfs will be generated
#  # setwd("/path/to/your/working/directory/")
#  
#  # And run:
#  artms_qualityControlEvidenceBasic(evidence_file = artms_data_ph_evidence,
#                                    keys_file = artms_data_ph_keys,
#                                    prot_exp =  "PH")

## ---- eval = FALSE---------------------------------------------------------
#  artms_writeConfigYamlFile(config_file_name = "config.yaml" )

## ---- eval = FALSE---------------------------------------------------------
#  artms_qualityControlEvidenceBasic(evidence_file = artms_data_ph_evidence,
#                                    keys_file = artms_data_ph_keys,
#                                    output_name = "qcPlots_evidence",
#                                    prot_exp = "PH")

## ---- eval = FALSE---------------------------------------------------------
#  artms_qualityControlEvidenceExtended(evidence_file = artms_data_ph_evidence,
#                                       keys_file = artms_data_ph_keys)

## ---- eval = FALSE---------------------------------------------------------
#  artms_qualityControlSummaryExtended(summary_file = "summary.txt",
#                                      keys_file = artms_data_ph_keys)

## ---- eval = FALSE---------------------------------------------------------
#  artms_quantification(
#    yaml_config_file = '/path/to/config/file/artms_ab_config.yaml')

## ---- eval = FALSE---------------------------------------------------------
#  artms_quantification(
#    yaml_config_file = '/path/to/config/file/artms_phglobal_config.yaml')

## ---- eval = FALSE---------------------------------------------------------
#  artms_proteinToSiteConversion(
#    evidence_file = "/path/to/the/evidence.txt",
#    ref_proteome_file = "/path/to/the/reference_proteome.fasta",
#    output_file = "/path/to/the/output/ph-sites-evidence.txt",
#    mod_type = "PH")

## ---- eval = FALSE---------------------------------------------------------
#  artms_proteinToSiteConversion(
#    evidence_file = "/path/to/the/evidence.txt",
#    ref_proteome_file = "/path/to/the/reference_proteome.fasta",
#    output_file = "/path/to/the/output/ub-sites-evidence.txt",
#    mod_type = "UB")

## ---- eval = FALSE---------------------------------------------------------
#  artms_quantification(
#    yaml_config_file = '/path/to/config/file/phsites_config.yaml')

## ---- echo = FALSE---------------------------------------------------------
#  artms_analysisQuantifications(log2fc_file = "ab-results.txt",
#                                modelqc_file = "ab-results_ModelQC.txt",
#                                specie = "human",
#                                output_dir = "AnalysisQuantifications")

