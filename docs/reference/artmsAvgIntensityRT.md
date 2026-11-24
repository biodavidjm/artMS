# Summarize average intensity and retention time per protein

Input an evidence file from MaxQuant and a file containing a list of
proteins of interest (optional). The function will summarize from the
evidence file and report back the average intensity, average retention
time, and the average caliberated retention time. If a list of proteins
is provided, then only those proteins will be summarized and returned.

## Usage

``` r
artmsAvgIntensityRT(
  evidence_file,
  protein_file = NULL,
  output_file = FALSE,
  verbose = TRUE
)
```

## Arguments

- evidence_file:

  (char) The filepath to the MaxQuant searched data (evidence) file (txt
  tab delimited file).

- protein_file:

  (char) The file path to a file or vector containing a list of proteins
  of interest.

- output_file:

  (char) The file name for the results (must have the extension `.txt`).
  If empty, then the results will be returned as an R object.

- verbose:

  (logical) `TRUE` (default) shows function messages

## Value

An R object with the results and a file with the results (if the
output_file argument is provided). It contains averages of Intensity,
Retention Time, Caliberated Retention Time

## Examples

``` r
ave_int <- artmsAvgIntensityRT(evidence_file = artms_data_ph_evidence)
#> >> READING IN FILES... 
#> --- Evidence file ready
#> >> COMPUTING AVERAGES... 
#> >> MERGING RESULTS... 
#> >> SUMMARIZATION COMPLETE!! 
```
