## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(
  eval=FALSE
)

## ---- eval = FALSE-------------------------------------------------------
#  # Evidence file
#  url_evidence <- '~/experiments/artms/ph/evidence.txt'
#  # url_evidence <- 'http://kroganlab.ucsf.edu/artms/ph/evidence.txt'
#  # evidence.df <- read.delim(url_evidence, stringsAsFactors = F)
#  
#  
#  # Keys file
#  url_keys <- "~/experiments/artms/ph/keys.txt"
#  # url_keys <- "http://kroganlab.ucsf.edu/artms/ph/keys.txt"
#  # keys.df <- read.delim(url_keys, stringsAsFactors = F)

## ---- eval = FALSE-------------------------------------------------------
#  artms_evidenceQC(evidence_file = url_evidence, keys_file = url_keys, prot_exp = 'ph')

## ---- eval = FALSE-------------------------------------------------------
#  artms_quantification(yaml_config_file = '/path/to/config/file/artms_ab_config.yaml')

## ---- eval = FALSE-------------------------------------------------------
#  artms_quantification(yaml_config_file = '/path/to/config/file/artms_phglobal_config.yaml')

## ---- eval = FALSE-------------------------------------------------------
#  artms_proteinToSiteConversion(evidence_file = "/path/to/the/evidence.txt", ref_proteome_file = "/path/to/the/reference_proteome.fasta", output_file = "/path/to/the/output/ph-sites-evidence.txt", mod_type = "PH")

## ---- eval = FALSE-------------------------------------------------------
#  artms_proteinToSiteConversion(evidence_file = "/path/to/the/evidence.txt", ref_proteome_file = "/path/to/the/reference_proteome.fasta", output_file = "/path/to/the/output/ub-sites-evidence.txt", mod_type = "UB")

## ---- eval = FALSE-------------------------------------------------------
#  artms_quantification(yaml_config_file = '/path/to/config/file/phsites_config.yaml')

