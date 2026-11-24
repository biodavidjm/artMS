# Check if a species is supported and available

Given a species name, it checkes whether is supported, and if supported,
check whether the annotation package is installed.

## Usage

``` r
artmsIsSpeciesSupported(species, verbose = TRUE)
```

## Arguments

- species:

  (char) The species name. Species currently supported as part of artMS:

  - HUMAN

  - MOUSE

  And the following species can be used as well, but the user needs to
  install the corresponding org.db package:

  - ANOPHELES (`install.packages(org.Ag.eg.db)`)

  - BOVINE (`install.packages(org.Bt.eg.db)`)

  - WORM (`install.packages(org.Ce.eg.db)`)

  - CANINE (`install.packages(org.Cf.eg.db)`)

  - FLY (`install.packages(org.Dm.eg.db)`)

  - ZEBRAFISH (`install.packages(org.Dr.eg.db)`)

  - CHICKEN (`install.packages(org.Gg.eg.db)`)

  - RHESUS (`install.packages(org.Mmu.eg.db)`)

  - CHIMP (`install.packages(org.Pt.eg.db)`)

  - RAT (`install.packages(org.Rn.eg.db)`)

  - YEAST (`install.packages(org.Sc.sgd.db)`)

  - PIG (`install.packages(org.Ss.eg.db)`)

  - XENOPUS (`install.packages(org.Xl.eg.db)`)

- verbose:

  (logical) `TRUE` (default) shows function messages

## Value

(string) Name of the package for the given species

## Examples

``` r
# Should return TRUE
artmsIsSpeciesSupported(species = "HUMAN")
#> [1] "org.Hs.eg.db"
artmsIsSpeciesSupported(species = "CHIMP")
#> ---(-) The package <org.Pt.eg.db> is not installed in your system.
#>            Just run: install.packages('org.Pt.eg.db') and try again
#> [1] FALSE
```
