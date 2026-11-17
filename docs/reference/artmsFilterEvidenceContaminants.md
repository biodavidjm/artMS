# Remove contaminants and empty proteins from the MaxQuant evidence file

Remove contaminants and erronously identified 'reverse' sequences by
MaxQuant, in addition to empty protein ids

## Usage

``` r
artmsFilterEvidenceContaminants(x, verbose = TRUE)
```

## Arguments

- x:

  (data.frame) of the Evidence file

- verbose:

  (logical) `TRUE` (default) shows function messages

## Value

(data.frame) without REV\_\_ and CON\_\_ Protein ids

## Examples

``` r
ef <- artmsFilterEvidenceContaminants(x = artms_data_ph_evidence)
#> -- Contaminants CON__|REV__ removed
```
