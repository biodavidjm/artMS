# Summarize the MSStats results and data quantification

Converts the MSStats results file to wide format (unique Protein ID and
columns are the comparisons), as well as adds BioReplicate information
about

- the Number of Unique Peptides,

- Spectral Counts

- Intensities for each protein. In cases where there are multiple values
  for a Protein-BioReplicate pair due to minute changes in sequence, the
  maximum value is taken for the pair. Any pairs without a value are
  assigned a value of NA.

## Usage

``` r
artmsMsstatsSummary(
  evidence_file,
  prot_group_file,
  keys_file,
  results_file,
  return_df = FALSE,
  verbose = TRUE
)
```

## Arguments

- evidence_file:

  (char or data.frame) The filepath to the MaxQuant searched data
  (evidence) file (txt tab delimited file). Only works for the newer
  versions of the evidence file.

- prot_group_file:

  (char) The filepath to the MaxQuant `proteinGroups.txt` file (txt tab
  delimited file) or data.frame

- keys_file:

  (char) The filepath to the keys file used with MSStats (txt tab
  delimited file).

- results_file:

  (char) The filepath to the MSStats results file in t he default long
  format (txt tab delimited file or data.frame).

- return_df:

  (data.frame) Whether or not to return the results to the R environment
  upon completion. This is useful if this is being used in an R pipeline
  and you want to feed the results directly into the next stage of
  analysis via an R environment/terminal. Regardless, the results will
  be written to file. Default = FALSE

- verbose:

  (logical) `TRUE` (default) shows function messages

## Value

(data.frame or txt file) with the summary

## Examples

``` r
# Testing warning if files are not submitted
test <- artmsMsstatsSummary(evidence_file = NULL,
                      prot_group_file = NULL,
                      keys_file = NULL,
                      results_file = NULL)
#> >> GENERATING A GLOBAL SUMMARY
```
