# MaxQuant evidence file to SAINTexpress format

Converts the MaxQuant evidence file to the 3 required files by
SAINTexpress. One can choose to either use the `spectral counts` (use
`msspc`) or the `intensities` (use `msint`) for the analysis.

## Usage

``` r
artmsEvidenceToSaintExpress(
  evidence_file,
  keys_file,
  ref_proteome_file,
  quant_variable = c("msspc", "msint"),
  output_file,
  verbose = TRUE
)
```

## Arguments

- evidence_file:

  (char) The evidence file path and name

- keys_file:

  (char) Keys file with a SAINT column specifying test (`T`) and control
  (`C`) conditions

- ref_proteome_file:

  (char) Reference proteome path file name in fasta format

- quant_variable:

  (char) choose either

  - `msspc` (spectral counts, default) or

  - `msint` (MS Intensity)

- output_file:

  (char) Output file name (must have extension .txt)

- verbose:

  (logical) `TRUE` (default) shows function messages

## Value

The 3 required files by SAINTexpress:

- `interactions.txt`

- `preys.txt`

- `baits.txt`

## Examples

``` r
# Testing that the files cannot be empty
artmsEvidenceToSaintExpress(evidence_file = NULL, 
keys_file = NULL, ref_proteome_file = NULL)
#> >> CONVERTING TO SAINTexpress FORMAT 
#> [1] "The evidence_file, keys_file, and ref_proteome_file \n           must not be empty"
```
