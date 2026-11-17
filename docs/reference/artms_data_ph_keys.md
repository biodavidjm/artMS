# Keys File Example

the `artMS` keys file provides the details of the experimental design
for any given proteomics experiment.

This particular example belongs to a PH experiment consisting of two
head and neck cancer cell lines ("Conditions" `"Cal33"` and `"HSC6"`),
with 2 biological replicates each (in this reduced version)

## Usage

``` r
artms_data_ph_keys
```

## Format

Tab delimited file with the following columns:

- Raw.file:

  Raw file processed. Each one should be a unique biological (or
  technical) replicate

- IsotopeLabelType:

  Type of labeling. `L` is used for label free experiments

- Condition:

  Label for conditions. VERY IMPORTANT: Only alpha-numeric characters
  and `underscore (_)` are allowed

- BioReplicate:

  Label for the Biological replicates. VERY IMPORTANT: Use the same
  labeling for bioreplicate as the Condition, but adding a `dash (-)`
  corresponding to the number of biological replicate. For example, for
  `Condition` `"Cal"`, use `Cal-1`, `Cal-2`, `Cal-3`, etc for the
  bioreplicates

- Run:

  The MS run number
