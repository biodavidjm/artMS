
# MANUSCRIPT DATASETS

######################################################
# Case Study: PHOSPHORYLATION
######################################################

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Dataset 1: PHOSPHORYLATION DANIELLE
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setwd('~/artMS_Manuscript/CaseStudy_PH_Danielle/')

# PH GLOBAL
artmsQuantification("phglobal/phglobal_config.yaml")

setwd("phglobal/")

artmsAnalysisQuantifications(log2fc_file = "phglobal-results.txt", 
                             modelqc_file = "phglobal-results_ModelQC.txt", 
                             species = "human",
                             output_dir = "analysisPhglobalTest")

setwd("analysisPhglobal_adjpvalue/")

artmsPhosfateOutput(inputFile = "phglobal-results-log2fc-long.txt", output_dir = "phosfate")

# PH SITES
artmsProtein2SiteConversion(evidence_file = "evidence.txt", 
                            ref_proteome_file = "uniprot_canonical.fasta", 
                            column_name = "Leading razor protein",
                            output_file = "phsites-evidence.txt",
                            mod_type = "PH")

artmsQuantification("phsites/phsites_config.yaml")

setwd("phsites/")

artmsAnalysisQuantifications(log2fc_file = "phsites-results.txt", 
                             modelqc_file = "phsites-results_ModelQC.txt", 
                             output_dir = "analysisQuant", 
                             isPtm = "ptmsites", 
                             species = "human")

artmsVolcanoPlot(mss_results = "phsites-results.txt")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Dataset 2: PHOSPHORYLATION PAPER OR8, OR9
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setwd("~/artMS_Manuscript/CaseStudy_PH_deGraaf_etal/")
evidence.ph <- read.delim("evidence.txt", stringsAsFactors = FALSE)
keys.ph <- read.delim("ph.keys.txt")
contrast.ph <- read.delim("ph.contrast.txt")

artmsProtein2SiteConversion(evidence_file = "evidence.txt", 
                            ref_proteome_file = "HUMAN_sprot-2012-09.fasta",
                            column_name = "Leading razor protein", 
                            output_file = "phsite-evidence.txt", 
                            mod_type = "PH")


# Create config files
# artmsWriteConfigYamlFile(config_file_name = "ph_global.yaml")
# artmsWriteConfigYamlFile(config_file_name = "ph_sites.yaml")

# PHGLOBAL

artmsQualityControlEvidenceBasic(evidence_file = "evidence.txt", 
                                 keys_file = "ph.option_1.keys.txt", 
                                 prot_exp = "PH", 
                                 output_name = "qc.phglobal.o1")

artmsQualityControlEvidenceExtended(evidence_file = "evidence.txt", 
                                    keys_file = "ph.option_1.keys.txt")


artmsQuantification(yaml_config_file = "phglobal.option_1.config.yaml")


setwd("phglobal_results_option1")

artmsAnalysisQuantifications(log2fc_file = "phglobal.results.txt", 
                             modelqc_file = "phglobal.results_ModelQC.txt", 
                             species = "human", 
                             output_dir = "analysisQuant")

setwd("../")

artmsQuantification(yaml_config_file = "ph.option_2.config.yaml")
setwd("results_phglobal_option2/")
artmsAnalysisQuantifications(log2fc_file = "phglobal.results.txt", 
                             modelqc_file = "phglobal.results_ModelQC.txt", 
                             species = "human", output_dir = "analysisQuant")


# checking normalized data
normalized <- read.delim("phglobal.results-mss-normalized.txt", stringsAsFactors = FALSE)

p1 <- ggplot2::ggplot(normalized, 
                      aes(x = GROUP_ORIGINAL, 
                          y = ABUNDANCE, 
                          fill = ABUNDANCE))
p1 <- p1 + geom_boxplot(aes(fill = GROUP_ORIGINAL))
p1 <- p1 + theme_linedraw()
p1 <-
  p1 + theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    legend.position = "none"
  )
p1 <- p1 + labs(x = "BIOREPLICATES")
p1 <- p1 + ggtitle("Relative Abundance BioReplicates")
print(p1)


artmsQuantification(yaml_config_file = "ph.OR8.config.yaml")


setwd("phglobal_OR8_results/")
artmsAnalysisQuantifications(log2fc_file = "phglobal.OR8.results.txt", 
                             modelqc_file = "phglobal.OR8.results_ModelQC.txt", 
                             output_dir = "analysisQuant", species = "human")





setwd("../")
artmsQuantification(yaml_config_file = "ph.OR9.config.yaml")

setwd("OR9_results_phglobal/")
artmsAnalysisQuantifications(log2fc_file = "phglobal.OR9.results.txt", 
                             modelqc_file = "phglobal.OR9.results_ModelQC.txt", 
                             output_dir = "analysisQuant", species = "human")




# PH SITES: 
artmsProtein2SiteConversion(evidence_file = "evidence.txt", 
                            ref_proteome_file = "HUMAN_sprot-2012-09.fasta", 
                            output_file = "phsite-evidence.txt", 
                            mod_type = "PH")

artmsQuantification(yaml_config_file = "ph.option_1.config.yaml")






##############################################################################
# CASE STUDY: CHANGES IN PROTEIN ABUNDANCE
##############################################################################

setwd("~/artMS_Manuscript/CaseStudy_ProteinAbundance/")

artmsWriteConfigYamlFile(config_file_name = "globalAbundance-config.yaml")

artmsQualityControlEvidenceBasic()

evidence_file = "ab-evidence.txt"
keys_file = "ab-keys.txt"
prot_exp =  "PH"
isSILAC = FALSE
plotINTDIST = TRUE
plotREPRO = TRUE
plotCORMAT = TRUE
plotINTMISC = TRUE
plotPTMSTATS = TRUE
printPDF = TRUE
verbose = TRUE

artmsQuantification(yaml_config_file = "globalAbundance-config.yaml")

setwd("results/")
artmsAnalysisQuantifications(log2fc_file = "ab-results.txt", 
                             modelqc_file = "ab-results_ModelQC.txt", 
                             species = "mouse", 
                             output_dir = "analysisQuant")


##############################################################################
# CASE STUDY: APMS
##############################################################################

# Regarding the controls, I prepared 4 replicates for each control control_HIVwt, 
# control_HIVwt_MG132 etc. However when I performed the final analysis 
# I combined the three baits (CUL5, ELOB and CBFB) for running the SAINT 
# analysis. Since I had 4 replicates of each control for each bait 
# (in total 12 replicates, I decided to use only 2 replicates that I 
#   prepared in parallel with each bait). Then I ran SAINT together 
# for all the three baits and had 6 controls for each condition.
# Please let me know if this makes sense. 

# The CUL5mock condition is not a control, only the ones that are labeled
# with control. CUL5mock condition was used to determine interactors of CUL5 
# in the absence of HIV infection. 
# So this is the way how I analyzed the data:
# 1. Use SAINT to determine specific interactions of CUL5 int he presence
# and absence of HIV infection
# 2. Using MSstats to determine differences in interactors between CUL5mock
# and CUL5wt.


setwd("~/artMS_Manuscript/CaseStudy_APMS_PXD009012/")

artmsQualityControlEvidenceBasic(evidence_file = "APMS-evidence.txt", 
                                 keys_file = "APMS-keys.txt", 
                                 prot_exp = "APMS")
artmsQualityControlEvidenceExtended(evidence_file = "APMS-evidence.txt", 
                                    keys_file = "APMS-keys.txt")





