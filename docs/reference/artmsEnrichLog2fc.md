# Enrichment of changes in protein abundance or PTMs

Enrichment analysis of the selected proteins

## Usage

``` r
artmsEnrichLog2fc(
  dataset,
  species,
  background,
  heatmaps = FALSE,
  output_name = "enrichment.txt",
  verbose = TRUE
)
```

## Arguments

- dataset:

  (data.frame) with a `Gene` and `Comparison or Label` (with the name of
  the comparisons specified in the contrast file) columns

- species:

  (char) Specie, only supported "human" or "mouse"

- background:

  (vector) Background genes for the enrichment analysis.

- heatmaps:

  (logical) if `TRUE` generates heatmaps (pdf), `FALSE` (default)
  otherwise.

- output_name:

  (char) Name of the annotation files, which will be used as well for
  the heatmaps (if `heatmaps` is selected) Default
  `output_name = "enrichment.txt"`

- verbose:

  (logical) `TRUE` (default) shows function messages

## Value

(data.frame) Results from the enrichment analysis using Gprofiler and
heatmaps (if selected)

## Examples

``` r
if (FALSE) { # \dontrun{
# The data must be annotated (Protein and Gene columns)
data_annotated <- artmsAnnotationUniprot(
                      x = artms_data_ph_msstats_results,
                      columnid = "Protein",
                      species = "human")
# And then the enrichment
enrich_set <- artmsEnrichLog2fc(
                   dataset = data_annotated,
                   species = "human",
                   background = unique(data_annotated$Gene))
} # }
```
