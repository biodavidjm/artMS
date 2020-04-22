#' artMS configuration template
#' 
#' @description The configuration file in `yaml` format contains 
#' the configuration details required to run `artmsQuantification()`, which
#' includes quality control functions
#' 
#' @format The configuration (`yaml`) file contains the following sections:
#' \describe{
#' \item{files}{ 
#' - `evidence` : /path/to/the/evidence.txt
#' - `keys` : /path/to/the/keys.txt
#' - `contrasts` : /path/to/the/contrast.txt
#' - `summary` : /path/to/the/summary.txt
#' - `output` : /path/to/the/output/results/results.txt}
#' 
#' \item{qc}{
#' - basic: 1 # 1 = yes; 0 = no
#' - extended: 1 # 1 = yes; 0 = no
#' - extendedSummary: 0 # 1 = yes; 0 = no
#' }
#' 
#' \item{data}{
#' - enabled : 1 # 1 = yes; 0 = no
#' - fractions: 
#'   - enabled : 0 # 1 for protein fractionation
#' - silac: 
#'   - enabled : 0 # 1 for SILAC experiments
#' - filters: 
#'      - enabled : 1
#' - contaminants : 1
#' - protein_groups : remove #remove, keep
#' - modifications : ab # PH, UB, AB, APMS
#' - sample_plots : 1 # correlation plots }
#' 
#' \item{msstats}{
#' -  enabled : 1
#' -  msstats_input : # blank if not previous msstats input file is available
#' -  profilePlots : none # before, after, before-after, none
#' -  normalization_method : equalizeMedians # globalStandards 
#' (include a reference protein(s) ), equalizeMedians, quantile, 0
#' -  normalization_reference :  #should be a value in the Protein column
#' -  summaryMethod : TMP # "TMP"(default) means Tukey's median polish, 
#' which is robust estimation method. "linear" uses linear mixed model. 
#' "logOfSum" conducts log2 (sum of intensities) per run.
#' -  censoredInt : NA  # Missing values are censored or at random. 'NA' 
#' (default) assumes that all 'NA's in 'Intensity' column are censored. 
#' '0' uses zero intensities as censored intensity. In this case, 
#' NA intensities are missing at random. The output from Skyline 
#' should use '0'. Null assumes that all NA intensites are randomly missing.
#' -  cutoffCensored : minFeature  # Cutoff value for censoring. only 
#' with censoredInt='NA' or '0'. Default is 'minFeature', which uses 
#' minimum value for each feature.'minFeatureNRun' uses the smallest between 
#' minimum value of corresponding feature and minimum value of corresponding
#'  run. 'minRun' uses minumum value for each run.
#' -  MBimpute : 1 # only for summaryMethod="TMP" and censoredInt='NA' or '0'.
#'  TRUE (default) imputes 'NA' or '0' (depending on censoredInt option) by 
#'  Accelated failure model. FALSE uses the values assigned by cutoffCensored.
#' -  feature_subset: all # all|highQuality  : highQuality seems to be buggy 
#' right now
#' }
#' 
#' \item{output_extras}{
#' - output_extras :
#'   - enabled : 1 # if 0, it wont do anything in this section
#' - annotate :  
#'   - enabled: 1 # 1|0 whether to annotate the proteins in the results or not
#' - species : HUMAN  # Supported species: HUMAN, MOUSE, ANOPHELES, ARABIDOPSIS, 
#' BOVINE, WORM, CANINE, FLY, ZEBRAFISH, ECOLI_STRAIN_K12, ECOLI_STRAIN_SAKAI, 
#' CHICKEN, RHESUS, MALARIA, CHIMP, RAT, YEAST, PIG, XENOPUS
#' - plots:
#'   - volcano: 1
#'   - heatmap: 1
#'   - LFC : -1.5 1.5 # Range of minimal log2fc
#'   - FDR : 0.05
#'   - heatmap_cluster_cols : 0
#'   - heatmap_display : log2FC # log2FC or pvalue}
#' }
"artms_config"


#' CORUM Protein Complexes database use for complex enrichment analysis
#' 
#' @description The list of protein complexes has been enriched with
#' mitochondria proteins from mouse, as described in this paper:
#' 
#' 2018 - Ruchi Masand, Esther Paulo, Dongmei Wu , Yangmeng Wang, 
#' Danielle L. Swaney, David Jimenez-Morales, Nevan J. Krogan, and Biao Wang
#' Proteome Imbalance of Mitochondrial Electron Transport Chain in Brown 
#' Adipocytes Leads to Metabolic Benefits. 
#' Cell Metab. 2018 Mar 06; 27(3):616-629.e4 
#' @details LAST CORUM DOWNLOAD DATE: 2017-08-01 
#' @format Tab delimited file.
#' \describe{
#' To find out more about the format and columns available at CORUM, 
#' please visit this 
#' [link](http://mips.helmholtz-muenchen.de/corum/)
#' }
"artms_data_corum_mito_database"

#' TB PATHOGEN: Mycobacterium tuberculosis 
#' (strain ATCC 35801 / TMC 107 / Erdman) UNIPROTS IDS
#'
#' @format A data.frame of Entry IDs
"artms_data_pathogen_TB"


#' LPN PATHOGEN: Legionella pneumophila subsp. pneumophila 
#' (strain Philadelphia 1 / ATCC 33152 / DSM 7513) UNIPROT IDS
#'
#' @format A data.frame of Entry IDs
"artms_data_pathogen_LPN"

#' Evidence file example
#'
#' @description Evidence file from a PH experiment consisting of two 
#' head and neck cancer cell lines ("Conditions" `"Cal33"` and `"HSC6"`). 
#' 
#' Unfortunately, the number of lines was reduced to 1/8 due 
#' to bioconductor limitations on data size, which means that this data is 
#' not very representative of a real evidence file. However, both the 
#' full evidence.txt and keys.txt file are available at:
#' http://kroganlab.ucsf.edu/artms/ph/evidence.txt
#' http://kroganlab.ucsf.edu/artms/ph/keys.txt
#'
#' @format A data frame with all the columns available in an evidence file
#' generated with MaxQuant version 1.6.2.3
"artms_data_ph_evidence"


#' Keys File Example
#' 
#' @description the `artMS` keys file provides the details of the experimental 
#' design for any given proteomics experiment. 
#' 
#' This particular example belongs to a PH experiment consisting of two 
#' head and neck cancer cell lines ("Conditions" `"Cal33"` and `"HSC6"`), 
#' with 2 biological replicates each (in this reduced version)
#' 
#' @format Tab delimited file with the following columns:
#' \describe{
#' \item{Raw.file}{Raw file processed. Each one should be a unique
#' biological (or technical) replicate}
#' 
#' \item{IsotopeLabelType}{Type of labeling. `L` is used for label free 
#' experiments}
#' 
#' \item{Condition}{Label for conditions. VERY IMPORTANT: Only alpha-numeric 
#' characters and `underscore (_)` are allowed}
#' 
#' \item{BioReplicate}{Label for the Biological replicates. VERY IMPORTANT:
#' Use the same labeling for bioreplicate as the Condition, but adding a 
#' `dash (-)` corresponding to the number of biological replicate. 
#' For example, for `Condition` `"Cal"`, use `Cal-1`, `Cal-2`, `Cal-3`, etc 
#' for the bioreplicates}
#' 
#' \item{Run}{The MS run number}
#' }
"artms_data_ph_keys"

#' Contrast example for the PH dataset
#'
#' @description Contrast file with the relative quantification to be performed
#' for the two conditions available in the example dataset: "Cal33-HSC6". 
#' See vignette for more details on how to prepare the contrast file.
#'
#' @format list with one comparison: "Cal33-HSC6"
"artms_data_ph_contrast"


#' MSstats results file example
#'
#' @description Relative quantification results obtained running MSstats
#' on a PH datasets (global analysis). Changes in protein phosphorylation
#' were quantified between two conditions
#'
#' @format A data frame resulting from running the lastest version of MSstats
"artms_data_ph_msstats_results"

#' Random data set
#'
#' Dataset randomly generated for testing purposes
#'
#' @format A data frame with 100 rows and 10 variables:
#' \describe{
#' Dataset generated using this code
#' 
#' `data.frame(replicate(10,sample(0:1,100,rep=TRUE)))`
#' }
"artms_data_randomDF"



#' artMS configuration for the available PH dataset
#' 
#' @description The configuration file with default options to run the
#' available PH dataset with `artmsQuantification()``
#' 
#' @format The configuration (`yaml`) file contains the following sections:
#' \describe{
#' \item{files}{ 
#' - `evidence` : artms_data_ph_evidence
#' - `keys` : artms_data_ph_keys
#' - `contrasts` : artms_data_ph_contrast
#' - `summary` : 
#' - `output` : "results.txt"
#' }
#' 
#' \item{qc}{
#' - basic: 0
#' - extended: 0
#' - extendedSummary: 0 =
#' }
#' 
#' \item{data}{
#' - enabled : 1 
#' - fractions: 
#'   - enabled : 0 
#' - silac: 
#'   - enabled : 0 
#' - filters: 
#'      - enabled : 1
#' - contaminants : 1
#' - protein_groups : remove 
#' - modifications : PH 
#' - sample_plots : 1 
#' }
#' 
#' \item{msstats}{
#' -  enabled : 1
#' -  msstats_input : # blank if not previous msstats input file is available
#' -  profilePlots : none # before, after, before-after, none
#' -  normalization_method : equalizeMedians 
#' -  normalization_reference :  #should be a value in the Protein column
#' -  summaryMethod : TMP 
#' -  censoredInt : NA
#' -  cutoffCensored : minFeature 
#' -  MBimpute : 1 
#' -  feature_subset: all 
#' }
#' 
#' \item{output_extras}{
#' - output_extras :
#'   - enabled : 1 
#' - annotate :  
#'   - enabled: 1 
#' - species : HUMAN  
#' - plots:
#'   - volcano: 1
#'   - heatmap: 1
#'   - LFC : -1 1 
#'   - FDR : 0.05
#'   - heatmap_cluster_cols : 0
#'   - heatmap_display : log2FC}
#' }
"artms_data_ph_config"


