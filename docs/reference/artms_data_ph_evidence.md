# Evidence file example

Evidence file from a PH experiment consisting of two head and neck
cancer cell lines ("Conditions" `"Cal33"` and `"HSC6"`).

Unfortunately, the number of lines was reduced to 1/20 due to
bioconductor limitations on data size, but it should be enough to test
the qc and quantification functions. The number of total columns from
the original evidence file was also reduced to 36 (out of the original
66 columns). Check `colnames(artms_data_ph_evidence)` for details

## Usage

``` r
artms_data_ph_evidence
```

## Format

A data frame with all the columns available in an evidence file
generated with MaxQuant version 1.6.2.3
