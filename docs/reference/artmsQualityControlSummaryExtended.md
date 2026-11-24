# Quality Control of the MaxQuant summary.txt file

Performs quality control based on the information available in the
MaxQuant summary.txt file.

## Usage

``` r
artmsQualityControlSummaryExtended(
  summary_file,
  keys_file,
  output_dir = "qc_summary",
  output_name = "qcExtended_summary",
  isFractions = FALSE,
  plotMS1SCANS = TRUE,
  plotMS2 = TRUE,
  plotMSMS = TRUE,
  plotISOTOPE = TRUE,
  printPDF = TRUE,
  verbose = TRUE
)
```

## Arguments

- summary_file:

  (char or data.frame) The evidence file path and name, or data.frame

- keys_file:

  (char or data.frame) The keys file path and name or data.frame

- output_dir:

  (char) Name for the folder to output the results plots. Default is
  "qc_summary".

- output_name:

  (char) prefix output name (no extension). Default:
  "qcExtended_summary"

- isFractions:

  (logical) `TRUE` if it is a 2D experiment (fractions). Default:
  `FALSE`

- plotMS1SCANS:

  (logical) `TRUE` generates MS1 scan counts plot: Page 1 shows the
  number of MS1 scans in each BioReplicate. If replicates are present,
  Page 2 shows the mean number of MS1 scans per condition with error bar
  showing the standard error of the mean. If isFractions `TRUE`, each
  fraction is a stack on the individual bar graphs.

- plotMS2:

  (logical) `TRUE` generates MS2 scan counts plot: Page 1 shows the
  number of MSs scans in each BioReplicate. If replicates are present,
  Page 2 shows the mean number of MS1 scans per condition with error bar
  showing the standard error of the mean. If isFractions `TRUE`, each
  fraction is a stack on the individual bar graphs.

- plotMSMS:

  (logical) `TRUE` generates MS2 identification rate (%) plot: Page 1
  shows the fraction of MS2 scans confidently identified in each
  BioReplicate. If replicates are present, Page 2 shows the mean rate of
  MS2 scans confidently identified per condition with error bar showing
  the standard error of the mean. If isFractions `TRUE`, each fraction
  is a stack on the individual bar graphs.

- plotISOTOPE:

  (logical) `TRUE` generates Isotope Pattern counts plot: Page 1 shows
  the number of Isotope Patterns with charge greater than 1 in each
  BioReplicate. If replicates are present, Page 2 shows the mean number
  of Isotope Patterns with charge greater than 1 per condition with
  error bar showing the standard error of the mean. If isFractions
  `TRUE`, each fraction is a stack on the individual bar graphs.

- printPDF:

  If `TRUE` (default) prints out the pdfs. Warning: plot objects are not
  returned due to the large number of them.

- verbose:

  (logical) `TRUE` (default) shows function messages

## Value

A number of plots from the summary file

## Examples

``` r
# Testing warning if files are not submitted
test <- artmsQualityControlSummaryExtended(summary_file = NULL,
keys_file = NULL)
#> ---------------------------------------------
#> artMS: EXTENDED QUALITY CONTROL (-summary.txt based)
#> ---------------------------------------------
```
