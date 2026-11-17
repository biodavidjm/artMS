# Relative quantification using MSstats

Relative quantification using MSstats including:

- plots

- quantifications (log2fc, pvalues, etc)

- normalized abundance values

## Usage

``` r
artmsQuantification(
  yaml_config_file,
  data_object = FALSE,
  printPDF = TRUE,
  printTables = TRUE,
  display_msstats = FALSE,
  return_results_object = FALSE,
  verbose = TRUE
)
```

## Arguments

- yaml_config_file:

  (char, required) The yaml file name and location

- data_object:

  (logical) flag to indicate whether the configuration file is a string
  to a file that should be opened or config object (yaml). Default is
  `FALSE`. Choose `TRUE` if `yaml_config_file` is a yaml object

- printPDF:

  (logical) if `TRUE` (default), prints out pdf

- printTables:

  (logical) `TRUE` (default) print results tables

- display_msstats:

  (logical) if `TRUE`, prints MSstats outputs (default is `FALSE`)

- return_results_object:

  (logical) Default is `FALSE`. If `TRUE`, it returns a list of data
  frames with MSstats results, including:

  - `comparisonResult`: comparison results

  - `ModelQC`

  - `FittedModel`: fit model details

  - `power`: power calculations

  - `sample_size`: sample size estimations

- verbose:

  (logical) `TRUE` (default) shows function messages

## Value

The relative quantification of the conditions and comparisons specified
in the keys/contrast file resulting from running MSstats, in addition to
quality control plots (if selected)

## Examples

``` r
# Recommended
# artmsQuantification(yaml_config_file = "your-config-file.yaml")

# Example to test this function using the example dataset available in artMS
# Step 1: Add evidence, keys, and contrast to configuration object
artms_data_ph_config$files$evidence <- artms_data_ph_evidence
artms_data_ph_config$files$keys <- artms_data_ph_keys
artms_data_ph_config$files$contrasts <- artms_data_ph_contrast

# Step 2: Run the quantification step
quant_results <- artmsQuantification(yaml_config_file = artms_data_ph_config, 
                                     data_object = TRUE, 
                                     display_msstats = FALSE,
                                     printPDF = FALSE,
                                     printTables = FALSE)
#> --------------------------------------------
#> artMS: Relative Quantification using MSstats
#> --------------------------------------------
#> >> Reading the configuration file
#> >> LOADING DATA 
#> >> MERGING FILES 
#> >> CONVERT Intensity values < 1 to NA
#> >> FILTERING 
#> -- Contaminants CON__|REV__ removed
#> -- Removing protein groups
#> -- Use <Leading.razor.protein> as Protein ID
#> >> CONVERTING THE DATA TO MSSTATS FORMAT 
#> -- Selecting Sequence Type: MaxQuant 'Modified.sequence' column
#>  (+) <Fraction> column added (with value 1, MSstats requirement)
#> -- Adding NA values for missing values (required by MSstats) 
#> >> RUNNING MSstats (it usually takes a 'long' time: please, be patient)
#>  (MSstats messages are turned off. Select <display_msstats = TRUE> to activate MSstats outputs)
#> >> MSstats done
# Check the list of data frames "quant_results". Nothing should be printed out. 
```
