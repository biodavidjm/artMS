# Individual Normalized abundance dot plots for every protein

Protein abundance dot plots for each unique uniprot id. It can take a
long time

## Usage

``` r
artmsDataPlots(input_file, output_file, verbose = TRUE)
```

## Arguments

- input_file:

  (char) File path and name to the `-normalized.txt` output file from
  MSstats

- output_file:

  (char) Output file (path) name (add the `.pdf` extension)

- verbose:

  (logical) `TRUE` (default) shows function messages

## Value

(pdf) file with each individual protein abundance plot for each
conditions

## Examples

``` r
if (FALSE) { # \dontrun{
artmsDataPlots(input_file = "results/ab-results-mss-normalized.txt",
               output_file = "results/ab-results-mss-normalized.pdf")
} # }
```
