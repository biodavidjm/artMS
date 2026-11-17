# MaxQuant evidence file to SAINTq format

Converts the MaxQuant evidence file to the required files by SAINTq. The
user can choose to use either peptides with `spectral counts` (use
`msspc`) or the all the peptides (use `all`) for the analysis. The
quantitative can be also chosen (either MS Intensity or Spectral Counts)

## Usage

``` r
artmsEvidenceToSAINTq(
  evidence_file,
  keys_file,
  output_dir = "artms_saintq",
  sc_option = c("all", "msspc"),
  fractions = FALSE,
  quant_variable = c("msint", "msspc"),
  verbose = TRUE
)
```

## Arguments

- evidence_file:

  (char or data.frame) The evidence file path and name, or data.frame

- keys_file:

  (char) Keys file with a SAINT column specifying test (`T`) and control
  (`C`) conditions

- output_dir:

  (char) New directory to create and save files. Default is current
  directory (recommended to provide a new folder name).

- sc_option:

  (char). Filter peptides with spectral counts only. Two options:

  - `msspc`: use only peptides with spectral_counts

  - `all` (default): all peptides detected (including the one resulting
    from the MaxQuant 'Match between run' algorithm)

- fractions:

  (logical) `TRUE` for 2D proteomics (fractions). Default: `FALSE`

- quant_variable:

  (char) Select the quantitative variable. Two options available:

  - `msint`: MS Intensity (default)

  - `msspc`: MS.MS.count (Spectral Counts)

- verbose:

  (logical) `TRUE` (default) shows function messages

## Value

The input files requires to run SAINTq

## Details

After running the script, the new specified folder should contain the
folling files:

- saintq-config-peptides

- saintq-config-proteins

- saintq_input_peptides.txt

- saintq_input_proteins.txt

Then `cd` into the new folder and run either of the following two
options (assuming that `saintq` is installed in your linux/unix/mac os x
system):

`> saintq config-saintq-peptides`

or

`> saintq config-saintq-proteins`

## Examples

``` r
# Testing that the files cannot be empty
artmsEvidenceToSAINTq   (evidence_file = NULL, 
                                   keys_file = NULL, 
                                   output_dir = NULL)
#> >> GENERATING A SAINTq INPUT FILE 
#> >> CHECKING THE keys FILE FIRST 
#> [1] "The evidence_file, keys_file and output_dir must not be NULL"
```
