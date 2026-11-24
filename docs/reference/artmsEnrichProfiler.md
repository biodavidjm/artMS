# Enrichment analysis using g:Profiler (gprofiler2)

This function simplifies the enrichment analysis performed by the
excellent tool g:Profiler, using the `gost()` function of the
`gprofiler2` package (which replaces the deprecated `gProfileR`
package). An internet connection is required.

## Usage

``` r
artmsEnrichProfiler(
  x,
  categorySource = c("GO"),
  species,
  background = NA,
  verbose = TRUE
)
```

## Arguments

- x:

  (list, data.frame) List of protein ids. It can be anything: either a
  list of ids, or you could also send a data.frame and it will find the
  columns with the IDs. Is not cool? Multiple list can be also sent
  simultaneously, as for example running:
  `tmp <- split(enrichment$Gene, enrichment$cl_number, drop= TRUE)`

- categorySource:

  (vector) Resources providing the terms on which the enrichment will be
  performed. The resources supported by g:Profiler are:

  - GO (GO:BP, GO:MF, GO:CC): Gene Ontology. `GO` selects the three of
    them

  - KEGG: Biological pathways

  - REAC: Biological pathways (Reactome)

  - WP: Biological pathways (WikiPathways)

  - TF: Regulatory motifs in DNA (TRANSFAC TFBS)

  - MIRNA: Regulatory motifs in DNA (miRTarBase microRNAs)

  - CORUM: protein complexes database

  - HP: Human Phenotype Ontology

  - HPA: Protein databases (Human Protein Atlas) `OMIM` is no longer
    available in g:Profiler and it is ignored. Electronic GO annotations
    (IEA) are excluded.

- species:

  (char) Specie code: Organism names are constructed by concatenating
  the first letter of the name and the family name. Example: human -
  ’hsapiens’, mouse - ’mmusculus’. Check g:Profiler to find out more
  about supported species.

- background:

  (vector) gene list to use as background for the enrichment analysis.
  Default: `NA` (no custom background: all known genes are used)

- verbose:

  (logical) `TRUE` (default) shows function messages

## Value

(data.frame) The significant enrichment results (zero rows if there are
none)

## Details

This function calls
[`gprofiler2::gost()`](https://rdrr.io/pkg/gprofiler2/man/gost.html)
with these arguments:

- ordered_query = FALSE

- significant = TRUE

- exclude_iea = TRUE

- measure_underrepresentation = FALSE

- evcodes = TRUE (needed to report the `intersection` column)

- user_threshold = 0.05

- correction_method = "g_SCS" (named "analytical" in gProfileR)

- domain_scope = "known", or "custom" when a `background` is provided

The column names of the results are kept as in previous versions of
artMS (based on `gProfileR`): `query.number`, `significant`, `p.value`,
`term.size`, `query.size`, `overlap.size`, `recall`, `precision`,
`term.id`, `domain`, `term.name`, `intersection`. Note that `domain` now
uses the g:Profiler source names (e.g. `GO:BP`, `KEGG`, `REAC`).

## Examples

``` r
if (FALSE) { # \dontrun{
# annotate the MSstats results to get the Gene name
data_annotated <- artmsAnnotationUniprot(
                                     x = artms_data_ph_msstats_results,
                                     columnid = "Protein",
                                     species = "human")

# Filter the list of genes with a log2fc > 2
filtered_data <- 
unique(data_annotated$Gene[which(data_annotated$log2FC > 2)])

# And perform enrichment analysis
data_annotated_enrich <- artmsEnrichProfiler(
                                   x = filtered_data,
                                   categorySource = c('KEGG'),
                                   species = "hsapiens",
                                   background = unique(data_annotated$Gene))
} # }
```
