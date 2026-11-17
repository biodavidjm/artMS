# Change a specific column name in a given data.frame

Making easier to change a column name in any data.frame

## Usage

``` r
artmsChangeColumnName(dataset, oldname, newname)
```

## Arguments

- dataset:

  (data.frame) with the column name you want to change

- oldname:

  (char) the old column name

- newname:

  (char) the new name for that column

## Value

(data.frame) with the new specified column name

## Examples

``` r
artms_data_ph_evidence <- artmsChangeColumnName(
                               dataset = artms_data_ph_evidence,
                               oldname = "Phospho..STY.",
                               newname = "PH_STY")
```
