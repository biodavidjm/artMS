#' artMS configuration template
#' 
#' @description The configuration file in `yaml` format contains 
#' the configuration details required to run `artms_quantification()`, which
#' includes quality control (`artms_evidenceQC()`) 
#' 
#' @format The configuration (`yaml`) file contains the following sections:
#' \describe{
#' \item{files}{ 
#' - `evidence` : /path/to/the/evidence.txt
#' - `keys` : /path/to/the/keys.txt
#' - `contrasts` : /path/to/the/contrast.txt
#' - `output` : /path/to/the/output/results/results.txt}
#' 
#' \item{qc}{
#' - `enabled` : 1 or 0
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
#' -  normalization_method : equalizeMedians # globalStandards (include a reference protein(s) ), equalizeMedians, quantile, 0
#' -  normalization_reference :  #should be a value in the Protein column
#' -  summaryMethod : TMP # "TMP"(default) means Tukey's median polish, which is robust estimation method. "linear" uses linear mixed model. "logOfSum" conducts log2 (sum of intensities) per run.
#' -  censoredInt : NA  # Missing values are censored or at random. 'NA' (default) assumes that all 'NA's in 'Intensity' column are censored. '0' uses zero intensities as censored intensity. In this case, NA intensities are missing at random. The output from Skyline should use '0'. Null assumes that all NA intensites are randomly missing.
#' -  cutoffCensored : minFeature  # Cutoff value for censoring. only with censoredInt='NA' or '0'. Default is 'minFeature', which uses minimum value for each feature.'minFeatureNRun' uses the smallest between minimum value of corresponding feature and minimum value of corresponding run. 'minRun' uses minumum value for each run.
#' -  MBimpute : 1 # only for summaryMethod="TMP" and censoredInt='NA' or '0'. TRUE (default) imputes 'NA' or '0' (depending on censoredInt option) by Accelated failure model. FALSE uses the values assigned by cutoffCensored.
#' -  feature_subset: all # all|highQuality  : highQuality seems to be buggy right now
#' }
#' 
#' \item{output_extras}{
#' - enabled : 0
#' - msstats_output : # keep empty if mstats is enabled. Or provide the mstats 'results.txt' file
#' - annotate : 0 # 1|0 whether to annotate the proteins in the results or not
#' - species : HUMAN  # can use multiple species, but separate with a "-" eg. HUMAN-MOUSE-HIV-...
#' - annotation_dir : /path/to/the/files
#' - comparisons : all # or any grep expression that returns a subset of the contrasts file
#' - LFC : -2 2
#' - FDR : 0.05
#' - heatmap : 1 
#' - heatmap_cluster_cols : 0
#' - heatmap_display : log2FC #or pvalue
#' - volcano : 1}
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

#' The artMS Keys File
#' 
#' @description the artMS Keys file provides the details of the experimental 
#' design  for any given proteomics experiment. 
#' 
#' To illustrate how the keys file works, the `ph_keys` data object
#' contains the information about a time series experiment from the 
#' [FluOMICS project](http://www.fluomics.org/). THP1 cells
#' were infected with Influenza and data collected at different time
#' points. in addition to the 
#' 
#' Data frame with the information about the experimental design.
#' It provides the metadata of the `ph_evidence` data object also available in 
#' this package (see `vignettes` to find out more)
#' 
#' Comparison of the phosphoproteome of two head and neck cancer cell lines, 
#' Cal33 and HSC6.
#' 
#' 4 replicates per cell line
#' 
#' @format Tab delimited file with the following columns:
#' \describe{
#' \item{Raw.file}{Raw file processed. Each one should be either a
#' biological (or technical) replicate}
#' 
#' \item{IsotopeLabelType}{Type of labeling. `L` is used for label free 
#' experiments}
#' 
#' \item{Condition}{Label for conditions. Only alpha-numeric characters and 
#' `underscore (_)` are allowed}
#' 
#' \item{BioReplicate}{Label for the Biological replicates. The consensus is
#' to use the same label as the Condition, but adding a `dash (-)` corresponding
#' to the number of biological replicate. For example, `Cal-1`, `Cal-2`, 
#' `Cal-3`}
#' 
#' \item{Run}{The MS run number}
#' }
"artms_data_keys_example"

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



