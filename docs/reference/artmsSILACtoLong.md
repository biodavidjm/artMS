# Convert the SILAC evidence file to MSstats format

Converting the evidence file from a SILAC search to a format compatible
with MSstats. It basically modifies the Raw.files adding the Heavy and
Light label

## Usage

``` r
artmsSILACtoLong(evidence_file, output = NULL, verbose = TRUE)
```

## Arguments

- evidence_file:

  (char) Text filepath to the evidence file

- output:

  (char) Text filepath of the output name. If NULL it does not write the
  output

- verbose:

  (logical) `TRUE` (default) shows function messages

## Value

(data.frame) with SILAC data processed for MSstats (and output file)

## Examples

``` r
if (FALSE) { # \dontrun{
evidence2silac <- artmsSILACtoLong(evidence_file = "silac.evicence.txt",
                                   output = "silac-evidence.txt")
} # }
```
