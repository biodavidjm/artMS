# artMS configuration for the available PH dataset

The configuration file with default options to run the available PH
dataset with \`artmsQuantification()“

## Usage

``` r
artms_data_ph_config
```

## Format

The configuration (`yaml`) file contains the following sections:

- files:

  - `evidence` : Empty. To test an example, run
    `artms_data_ph_config$files$evidence <- artms_data_ph_evidence`

  - `keys` : Empty To test an example datasets run
    `artms_data_ph_config$files$keys <- artms_data_ph_keys`

  - `contrasts` : Empty. To test the example datasets, run
    `artms_data_ph_config$files$contrasts <- artms_data_ph_contrast`

  - `summary` :

  - `output` : "results.txt"

- qc:

  - basic: 0

  - extended: 0

  - extendedSummary: 0

- data:

  - enabled : 1

  - silac:

    - enabled : 0

  - filters:

    - enabled : 1

  - contaminants : 1

  - protein_groups : remove

  - modifications : PH

  - sample_plots : 1

- msstats:

  - enabled : 1

  - msstats_input : \# blank if not previous msstats input file is
    available

  - profilePlots : none \# before, after, before-after, none

  - normalization_method : equalizeMedians

  - normalization_reference : \#should be a value in the Protein column

  - summaryMethod : TMP

  - censoredInt : NA

  - cutoffCensored : minFeature

  - MBimpute : 1

  - feature_subset: all

- output_extras:

  - output_extras :

    - enabled : 1

  - annotate :

    - enabled: 1

  - species : HUMAN

  - plots:

    - volcano: 1

    - heatmap: 1

    - LFC : -1 1

    - FDR : 0.05

    - heatmap_cluster_cols : 0

    - heatmap_display : log2FC
