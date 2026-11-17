# Write out a template file of the artMS configuration file (yaml)

Creates a template file of the artMS configuration file, which is
required to run `artmsQuantification`. Check
[`?artms_config`](http://artms.org/reference/artms_config.md) and the
vignettes to find out more about the details of the structure of the
file and how to fill it up

## Usage

``` r
artmsWriteConfigYamlFile(
  config_file_name = "artms_config_file.yaml",
  overwrite = FALSE,
  verbose = TRUE
)
```

## Arguments

- config_file_name:

  (char) The name for the configuration file. It must have a `.yaml`
  extension. If `NULL`, it returns the config as a yaml object

- overwrite:

  (logical) Default FALSE

- verbose:

  (logical) `TRUE` (default) shows function messages

## Value

A file (or yaml data object) of the artMS configuration file

## Examples

``` r
config_empty <- artmsWriteConfigYamlFile(config_file_name = NULL)
```
