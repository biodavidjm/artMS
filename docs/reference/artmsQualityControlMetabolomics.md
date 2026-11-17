# Quality Control analysis of the evidence-like metabolomics dataset

Quality Control analysis of the evidence-like metabolomics dataset

## Usage

``` r
artmsQualityControlMetabolomics(
  evidence_file,
  keys_file,
  met_exp = c("MV"),
  output_name = "qcPlots_metab",
  plotINTDIST = FALSE,
  plotCORMAT = TRUE,
  plotINTMISC = TRUE,
  printPDF = TRUE,
  verbose = TRUE
)
```

## Arguments

- evidence_file:

  (char or data.frame) The evidence file path and name, or data.frame

- keys_file:

  (char or data.frame) The keys file path and name or data.frame

- met_exp:

  (char) Metabolomics experiment. Only one option available (so far):

  - `MV`: Markview output

- output_name:

  (char) prefix output name (no extension). Default: "qcPlots_metab"

- plotINTDIST:

  if `TRUE` (default) plots both *Box-dot plot* and *Jitter plot* of
  biological replicates based on MS (raw) intensity values.

- plotCORMAT:

  if `TRUE` (default) generates up to 3 pdf files for technical
  replicates, biological replicates, and conditions. Each pdf file
  contains:

  - *Correlation matrix* for all the biological replicates using MS
    Intensity values,

  - *Clustering matrix* of the MS Intensities and correlation
    distribution

  - *histogram* of the distribution of correlations

- plotINTMISC:

  if `TRUE` (default) plots several pages, including bar plots of *Total
  Sum of Intensities in BioReplicates*, *Total Sum of Intensities in
  Conditions*, *Total Feature Counts in BioReplicates*, *Total Feature
  Counts in conditions* separated by categories (INT: has a intensity
  value NOINT: no intensity value ) *Box plots* of MS Intensity values
  per biological replicates and conditions; *bar plots* of total
  intensity by bioreplicates and conditions; Barplots of *total feature
  counts* by bioreplicates and conditions.

- printPDF:

  If `TRUE` (default) prints out the pdfs. Warning: plot objects are not
  returned due to the large number of them.

- verbose:

  (logical) `TRUE` (default) shows function messages

## Value

Quality control files and plots for metabolomics

## Examples

``` r
# Testing that input arguments cannot be null
artmsQualityControlMetabolomics(evidence_file = NULL,
                 keys_file = NULL,
                 met_exp = "MV")
#> [1] "Both <evidence_file> and <keys_file> cannot be NULL"
```
