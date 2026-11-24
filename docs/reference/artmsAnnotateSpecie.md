# Adding a column with the species name

Adding the species name to every protein. This makes more sense if there
are more than one species in the dataset, which must be specified in the
`pathogen` option. Influenza is a special case that it does not need to
be specified, as far as the proteins were originally annotated as
`INFLUENZAGENE_STRAIN` (strains covered `H1N1`, `H3N2`, `H5N1`), as for
example, `NS1_H1N1`

## Usage

``` r
artmsAnnotateSpecie(df, pathogen = "nopathogen", species)
```

## Arguments

- df:

  (data.frame) with a `Protein` column (of uniprot ids)

- pathogen:

  (char) Is there a pathogen in the dataset as well? if it does not,
  then use `pathogen = nopathogen` (default). Supported`tb`
  (Tuberculosis), `lpn` (Legionella)

- species:

  (char) Host organism (supported for now: `human` or `mouse`)

## Value

(data.frame) The same data.frame but with an extra column specifying the
species

## Examples

``` r
# Adding a new column with the main species of the data. Easy.
# But the main functionality is to add both the host-species and a pathogen,
# which is not illustrated in this example
data_with_specie <- artmsAnnotateSpecie(df = artms_data_ph_msstats_results,
                                         species = "human")
```
