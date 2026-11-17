# Outputs a heatmap of the MSStats results created using the log2fold changes

Heatmap of the Relative Quantifications (MSStats results)

## Usage

``` r
artmsPlotHeatmapQuant(
  input_file,
  output_file = "quantifications_heatmap.pdf",
  species,
  labels = "*",
  cluster_cols = FALSE,
  display = "log2FC",
  lfc_lower = -2,
  lfc_upper = 2,
  whatPvalue = "adj.pvalue",
  FDR = 0.05,
  verbose = TRUE
)
```

## Arguments

- input_file:

  (char) MSstats `results.txt` file and location (or data.frame of
  resuts)

- output_file:

  (char) Output file name (pdf format) and location.
  Default:"quantifications_heatmap.pdf"

- species:

  (char). Specie name to be able to add the Gene name. To find out more
  about the supported species check
  [`?artmsMapUniprot2Entrez`](http://artms.org/reference/artmsMapUniprot2Entrez.md)

- labels:

  (vector) of uniprot ids if only specific labes would like to be
  plotted. Default: all labels

- cluster_cols:

  (boolean) `True` or `False` to cluster columns. Default: FALSE

- display:

  Metric to be displayed. Options:

  - `log2fc` (default)

  - `adj.pvalue`

  - `pvalue`

- lfc_lower:

  (int) Lower limit for the log2fc. Default: -2

- lfc_upper:

  (int) Upper limit for the log2fc. Default: +2

- whatPvalue:

  (char) `pvalue` or `adj.pvalue` (default)

- FDR:

  (int) Upper limit false discovery rate (or pvalue). Default: 0.05

- verbose:

  (logical) `TRUE` (default) shows function messages

## Value

(pdf or ggplot2 object) heatmap of the MSStats results using the
selected metric

## Examples

``` r
# Unfortunately, the example does not contain any significant hits
# Use for illustration purposes
artmsPlotHeatmapQuant(input_file = artms_data_ph_msstats_results,
                       species = "human",
                       output_file = NULL,
                       whatPvalue = "pvalue",
                       lfc_lower = -1,
                       lfc_upper = 1)
#> >> TOTAL NUMBER OF SELECTED HITS FOR PLOTS WITH LFC BETWEEN -1 AND 1 AT 0.05 FDR:984
```
