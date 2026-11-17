# Generate ph-site specific detailed file

Generate extended detailed ph-site file, where every line is a ph site
instead of a peptide. Therefore, if one peptide has multiple ph sites it
will be breaking down in each of the sites. This file will help generate
input files for tools as [Phosfate](http://phosfate.com/) or
[PHOTON](https://github.com/jdrudolph/photon)

## Usage

``` r
artmsGeneratePhSiteExtended(
  df,
  pathogen = "nopathogen",
  species,
  ptmType,
  output_name
)
```

## Arguments

- df:

  (data.frame) of log2fc and imputed values

- pathogen:

  (char) Is there a pathogen in the dataset as well? Available pathogens
  are `tb` (Tuberculosis), `lpn` (Legionella). If it is not, then use
  `nopathogen` (default).

- species:

  (char) Main organism (supported for now: `human` or `mouse`)

- ptmType:

  (char) It must be a ptm-site quantification dataset. Either: yes:
  `ptmsites` (for site specific analysis), or `ptmph` (Jeff's script
  output evidence file).

- output_name:

  (char) A output file name (extension `.txt` required)

## Value

(data.frame) extended version of the ph-site

## Examples

``` r
if (FALSE) { # \dontrun{
artmsGeneratePhSiteExtended(df = dfobject, 
                             species = "mouse", 
                             ptmType = "ptmsites",
                             output_name = log2fc_file)
} # }
```
