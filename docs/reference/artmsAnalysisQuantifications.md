# Analysis of the Relative Quantifications

Analysis of relative quantifications, including:

- Annotations

- Summary files in different format (xls, txt) and shapes (long, wide)

- Numerous summary plots

- Enrichment analysis using Gprofiler

- PCA of quantifications

- Clustering analysis

- Basic imputation of missing values

To run this function, the following packages must be installed on your
system:

- From bioconductor:
  `BiocManager::install(c("ComplexHeatmap", "org.Mm.eg.db"))`

- From CRAN:
  `install.packages(c("factoextra", "FactoMineR", "gprofiler2", "PerformanceAnalytics"))`

## Usage

``` r
artmsAnalysisQuantifications(
  log2fc_file,
  modelqc_file,
  species,
  output_dir = "analysis_quant",
  outliers = c("keep", "iqr", "std"),
  enrich = TRUE,
  l2fc_thres = 1,
  choosePvalue = c("adjpvalue", "pvalue"),
  isBackground = "nobackground",
  isPtm = "global",
  mnbr = 2,
  pathogen = "nopathogen",
  plotPvaluesLog2fcDist = TRUE,
  plotAbundanceStats = TRUE,
  plotReproAbundance = TRUE,
  plotCorrConditions = TRUE,
  plotCorrQuant = TRUE,
  plotPCAabundance = TRUE,
  plotFinalDistributions = TRUE,
  plotPropImputation = TRUE,
  plotHeatmapsChanges = TRUE,
  plotTotalQuant = TRUE,
  plotClusteringAnalysis = TRUE,
  data_object = FALSE,
  printPDF = TRUE,
  verbose = TRUE
)
```

## Arguments

- log2fc_file:

  (char) MSstats results file location

- modelqc_file:

  (char) MSstats modelqc file location

- species:

  (char) Select one species. Species currently supported for a full
  analysis (including enrichment analysis):

  - HUMAN

  - MOUSE To find out species supported only for annotation check
    `?artmsIsSpeciesSupported()`

- output_dir:

  (char) Name for the folder to output the results from the function.
  Default is current directory (recommended to provide a new folder
  name).

- outliers:

  (char) It allows to keep or remove outliers. Options:

  - `keep` (default): it keeps outliers 'keep', 'iqr', 'std'

  - `iqr` (recommended): remove outliers +/- 6 x Interquartile Range
    (IQR)

  - `std` : 6 x standard deviation

- enrich:

  (logical) Performed enrichment analysis using GprofileR? Only
  available for species HUMAN and MOUSE. `TRUE` (default if "human" or
  "mouse" are the species) or `FALSE`

- l2fc_thres:

  (int) log2fc cutoff for enrichment analysis (default,
  `l2fc_thres = 1.5`)

- choosePvalue:

  (char) specify whether `pvalue` or `adjpvalue` should use for the
  analysis. The default option is `adjpvalue` (multiple testing
  correction). But if the number of biological replicates for a given
  experiment is too low (for example n = 2), then
  `choosePvalue = pvalue` is recommended.

- isBackground:

  (char) background of gene names for enrichment analysis.
  `nobackground` (default) will use the total number of genes detected.
  Alternatively provided the file path name to the background gene list.

- isPtm:

  (char) Is a ptm-site quantification?

  - `global` (default),

  - `ptmsites` (for site specific analysis),

  - `ptmph` (Jeff Johnson script output evidence file)

- mnbr:

  (int) **PARAMETER FOR NAIVE IMPUTATION**: "minimal number of
  biological replicates" for "naive imputation" and filtering. Default:
  `mnbr = 2`. Details: Intensity values for proteins/PTMs that are
  completely missed in one of the two conditions compared ("condition
  A"), but are found in at least 2 biological replicates (`mnbr = 2`) of
  the other "condition B", are imputed (values artificially assigned)
  and the log2FC values calculated. The goal is to keep those
  proteins/PTMs that are consistently found in one of the two conditions
  (in this case "condition B") and facilitate the inclusion in
  downstream analysis (if wished). The imputed intensity values are
  sampled from the lowest intensity values detected in the experiment,
  and (**WARNING**) the p-values are just randomly assigned between 0.05
  and 0.01 for illustration purposes (when generating a volcano plot
  with the output of `artmsAnalysisQuantifications`) or to include them
  when making a cutoff of p-value \< 0.05 for enrichment analysis
  **CAUTION**: `mnbr` would also add the constraint that any protein
  must be identified in at least `nmbr` biological replicates of the
  **same** condition or it will be filtered out. That is, if `mnbr = 2`,
  a protein found in two conditions but only in one biological replicate
  in each of them, it would be removed.

- pathogen:

  (char) Is there a pathogen in the dataset as well? if it does not,
  then use `pathogen = nopathogen` (default). Pathogens available: `tb`
  (Tuberculosis), `lpn` (Legionella)

- plotPvaluesLog2fcDist:

  (logical) If `TRUE` (default) plots pvalues and log2fc distributions

- plotAbundanceStats:

  (logical) If `TRUE` (default) plots stats graphs about abundance
  values

- plotReproAbundance:

  (logical) If `TRUE` plots reproducibility based on normalized
  abundance values

- plotCorrConditions:

  (logical) If `TRUE` plots correlation between the different conditions

- plotCorrQuant:

  (logical) if `TRUE` plots correlation between the available
  quantifications (comparisons)

- plotPCAabundance:

  (logical) if `TRUE` performs PCA analysis of conditions using
  normalized abundance values

- plotFinalDistributions:

  (logical) if `TRUE` plots distribution of both log2fc and pvalues

- plotPropImputation:

  (logical) if `TRUE` plots proportion of overall imputation

- plotHeatmapsChanges:

  (logical) if `TRUE` plots heatmaps of quantified changes (both all and
  significant only). Only if `printPDF` is also `TRUE`

- plotTotalQuant:

  (logical) if `TRUE` plots barplot of total number of quantifications
  per comparison

- plotClusteringAnalysis:

  (logical) if `TRUE` performs clustering analysis between quantified
  comparisons (more than 1 comparison required)

- data_object:

  (logical) flag to indicate whether the required files are data
  objects. Default is FALSE

- printPDF:

  If `TRUE` (default) prints out the pdfs. Warning: plot objects are not
  returned due to the large number of them.

- verbose:

  (logical) `TRUE` (default) shows function messages

## Value

(data.frame) summary of quantifications, including annotations,
enrichments, etc

## Examples

``` r
# Testing that the files cannot be empty
artmsAnalysisQuantifications(log2fc_file = NULL,
                              modelqc_file = NULL,
                              species = NULL,
                              output_dir = NULL)
#> ---------------------------------------------
#> artMS: ANALYSIS OF QUANTIFICATIONS
#> ---------------------------------------------
#> [1] "The evidence_file, modelqc_file, species and output_dir arguments cannot be NULL"
```
