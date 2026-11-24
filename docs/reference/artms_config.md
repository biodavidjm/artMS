# artMS configuration template

The configuration file in `yaml` format contains the configuration
details required to run
[`artmsQuantification()`](http://artms.org/reference/artmsQuantification.md),
which includes quality control functions

## Usage

``` r
artms_config
```

## Format

The configuration (`yaml`) file contains the following sections:

- files:

  - `evidence` : /path/to/the/evidence.txt

  - `keys` : /path/to/the/keys.txt

  - `contrasts` : /path/to/the/contrast.txt

  - `summary` : /path/to/the/summary.txt

  - `output` : /path/to/the/output/results/results.txt

- qc:

  - basic: 1 \# 1 = yes; 0 = no

  - extended: 1 \# 1 = yes; 0 = no

  - extendedSummary: 0 \# 1 = yes; 0 = no

- data:

  - enabled : 1 \# 1 = yes; 0 = no

  - silac:

    - enabled : 0 \# 1 for SILAC experiments

  - filters:

    - enabled : 1

  - contaminants : 1

  - protein_groups : remove \#remove, keep

  - modifications : ab \# PH, UB, AB, APMS

  - sample_plots : 1 \# correlation plots

- msstats:

  - enabled : 1

  - msstats_input : \# blank if not previous msstats input file is
    available

  - profilePlots : none \# before, after, before-after, none

  - normalization_method : equalizeMedians \# globalStandards (include a
    reference protein(s) ), equalizeMedians, quantile, 0

  - normalization_reference : \#should be a value in the Protein column

  - summaryMethod : TMP \# "TMP"(default) means Tukey's median polish,
    which is robust estimation method. "linear" uses linear mixed model.
    "logOfSum" conducts log2 (sum of intensities) per run.

  - censoredInt : NA \# Missing values are censored or at random. 'NA'
    (default) assumes that all 'NA's in 'Intensity' column are censored.
    '0' uses zero intensities as censored intensity. In this case, NA
    intensities are missing at random. The output from Skyline should
    use '0'. Null assumes that all NA intensites are randomly missing.

  - MBimpute : 1 \# only for summaryMethod="TMP" and censoredInt='NA' or
    '0'. TRUE (default) imputes 'NA' or '0' (depending on censoredInt
    option) by Accelated failure model. FALSE uses the values assigned
    by cutoffCensored

  - For all othe features, please check documentation for MSstats'
    dataProcess function

- output_extras:

  - output_extras :

    - enabled : 1 \# if 0, it wont do anything in this section

  - annotate :

    - enabled: 1 \# 1\|0 whether to annotate the proteins in the results
      or not

  - species : HUMAN \# Supported species: HUMAN, MOUSE, ANOPHELES,
    ARABIDOPSIS, BOVINE, WORM, CANINE, FLY, ZEBRAFISH, ECOLI_STRAIN_K12,
    ECOLI_STRAIN_SAKAI, CHICKEN, RHESUS, MALARIA, CHIMP, RAT, YEAST,
    PIG, XENOPUS

  - plots:

    - volcano: 1

    - heatmap: 1

    - LFC : -1.5 1.5 \# Range of minimal log2fc

    - FDR : 0.05

    - heatmap_cluster_cols : 0

    - heatmap_display : log2FC \# log2FC or pvalue
