# Reshape the MSstats results file from long to wide format

Converts the normal MSStats results.txt file into "wide" format where
each row represents a unique protein's results, and each column
represents the comparison made by MSStats. The fold change and p-value
of each comparison will be its own column.

## Usage

``` r
artmsResultsWide(
  results_msstats,
  output_file = NULL,
  select_pvalues = c("adjpvalue", "pvalue"),
  species,
  verbose = TRUE
)
```

## Arguments

- results_msstats:

  (char) Input file name and location (MSstats `results.txt` file)

- output_file:

  (char) Output file name and location (e.g. `results-wide.txt`). If
  `NULL` (default) returns an R object (data.frame)

- select_pvalues:

  (char) Either

  - `pvalue` or

  - `adjpvalue` (default)

- species:

  (char) Specie name for annotation purposes. Check
  [`?artmsMapUniprot2Entrez`](http://artms.org/reference/artmsMapUniprot2Entrez.md)
  to find out more about the supported species (e.g `species = "human"`)

- verbose:

  (logical) `TRUE` (default) shows function messages

## Value

(output file tab delimited) reshaped file with unique protein ids and as
many columns log2fc and adj.pvalues as comparisons available

## Examples

``` r
ph_results_wide <- artmsResultsWide(
                         results_msstats = artms_data_ph_msstats_results,
                         output_file = NULL,
                         species = "human")
#> >> RESHAPING MSSTATS RESULTS TO wide FORMAT 
```
