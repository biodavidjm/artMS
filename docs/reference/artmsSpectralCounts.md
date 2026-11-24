# Outputs the spectral counts from the MaxQuant evidence file.

Outputs the spectral counts from the MaxQuant evidence file.

## Usage

``` r
artmsSpectralCounts(
  evidence_file,
  keys_file,
  output_file = NULL,
  verbose = TRUE
)
```

## Arguments

- evidence_file:

  (char) Maxquant evidence file or data object

- keys_file:

  (char) Keys file with the experimental design or data object

- output_file:

  (char) Output file name (add `.txt` extension). If `NULL` (default) it
  returns a data.frame object

- verbose:

  (logical) `TRUE` (default) shows function messages

## Value

A txt file with biological replicates, protein id, and spectral count
columns

## Examples

``` r
summary_spectral_counts <- artmsSpectralCounts(
                                 evidence_file = artms_data_ph_evidence,
                                 keys_file = artms_data_ph_keys)
#> >> EXTRACTING SPECTRAL COUNTS FROM THE EVIDENCE FILE 
#> >> MERGING FILES 
```
