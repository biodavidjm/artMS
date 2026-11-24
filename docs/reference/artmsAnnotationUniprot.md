# Annotate table with Gene Symbol and Name based on Uniprot ID(s)

Annotate gene name and symbol based on uniprot ids. It will take the
column from your data.frame specified by the `columnid` argument, search
for the gene symbol, name, and entrez based on the species (`species`
argument) and merge the information back to the input data.frame

## Usage

``` r
artmsAnnotationUniprot(x, columnid, species, verbose = TRUE)
```

## Arguments

- x:

  (data.frame) to be annotated (or file path and name)

- columnid:

  (char) The column with the uniprotkb ids

- species:

  (char) The species name. Check
  [`?artmsMapUniprot2Entrez`](http://artms.org/reference/artmsMapUniprot2Entrez.md)
  to find out more about supported species.

- verbose:

  (logical) `TRUE` (default) shows function messages

## Value

(data.frame) with two new columns: `Gene` and `Protein.name`

## Examples

``` r
# This example adds annotations to the example evidence file included in
# artMS, based on the column 'Proteins'.

evidence_anno <- artmsAnnotationUniprot(x = artms_data_ph_evidence,
                                         columnid = 'Proteins',
                                         species = 'human')
```
