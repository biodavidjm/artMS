## ---- echo = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  eval=FALSE
)

## ---- eval = FALSE---------------------------------------------------------
#  suppressMessages(library(artMS))
#  artms_evidenceQCbasic(evidence_file = artms_data_ph_evidence,
#                   keys_file = artms_data_ph_keys,
#                   prot_exp = 'PH')

## ---- eval = FALSE---------------------------------------------------------
#  artms_quantification(yaml_config_file = '/path/to/config/file/artms_ab_config.yaml')

## ---- eval = FALSE---------------------------------------------------------
#  artms_quantification(yaml_config_file = '/path/to/config/file/artms_phglobal_config.yaml')

## ---- eval = FALSE---------------------------------------------------------
#  artms_proteinToSiteConversion(evidence_file = "/path/to/the/evidence.txt", ref_proteome_file = "/path/to/the/reference_proteome.fasta", output_file = "/path/to/the/output/ph-sites-evidence.txt", mod_type = "PH")

## ---- eval = FALSE---------------------------------------------------------
#  artms_proteinToSiteConversion(evidence_file = "/path/to/the/evidence.txt", ref_proteome_file = "/path/to/the/reference_proteome.fasta", output_file = "/path/to/the/output/ub-sites-evidence.txt", mod_type = "UB")

## ---- eval = FALSE---------------------------------------------------------
#  artms_quantification(yaml_config_file = '/path/to/config/file/phsites_config.yaml')

