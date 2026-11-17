# Map GENE SYMBOL, NAME, AND ENTREZID to a vector of Uniprot IDS

Map GENE SYMBOL, NAME, AND ENTREZID to a vector of Uniprot IDS

## Usage

``` r
artmsMapUniprot2Entrez(uniprotkb, species)
```

## Arguments

- uniprotkb:

  (vector) Vector of UniprotKB IDs

- species:

  (char) The species name. Species currently supported as part of artMS:
  check `?artmsIsSpeciesSupported()` to find out the list of supported
  species\`

## Value

(data.frame) with ENTREZID and GENENAMES mapped on UniprotKB ids

## Examples

``` r
# Load an example with human proteins
exampleID <- c("Q6P996", "B1N8M6")
artmsMapUniprot2Entrez(uniprotkb = exampleID, 
                       species = "HUMAN")
#>   UNIPROT  SYMBOL
#> 1  Q6P996  PDXDC1
#> 2  B1N8M6 KIR2DS1
#>                                                                                GENENAME
#> 1                                 pyridoxal dependent decarboxylase domain containing 1
#> 2 killer cell immunoglobulin like receptor, two Ig domains and short cytoplasmic tail 1
#>   ENTREZID
#> 1    23042
#> 2     3806
```
