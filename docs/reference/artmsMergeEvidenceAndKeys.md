# Merge evidence.txt (or summary.txt) with keys.txt files

Merge the evidence and keys files on the given columns

## Usage

``` r
artmsMergeEvidenceAndKeys(
  x,
  keys,
  by = c("RawFile"),
  isSummary = FALSE,
  verbose = TRUE
)
```

## Arguments

- x:

  (data.frame or char) The evidence data, either as data.frame or the
  file name (and path). It also works for the summary.txt file

- keys:

  The keys data, either as a data.frame or file name (and path)

- by:

  (vector) specifying the columns use to merge the evidence and keys.
  Default: `by=c('RawFile')`

- isSummary:

  (logical) TRUE or FALSE (default)

- verbose:

  (logical) `TRUE` (default) shows function messages

## Value

(data.frame) with the evidence and keys merged

## Examples

``` r
evidenceKeys <- artmsMergeEvidenceAndKeys(x = artms_data_ph_evidence,
                                           keys = artms_data_ph_keys)
#> >> MERGING FILES 
```
