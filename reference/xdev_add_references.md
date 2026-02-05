# xdev_add_references

Add resource references to a
[dataset](https://mothur.org/strollur/reference/dataset.md) object

## Usage

``` r
xdev_add_references(
  data,
  table,
  reference_name = "reference_names",
  reference_version = "reference_versions",
  reference_usage = "reference_usages",
  reference_note = "reference_notes",
  reference_url = "reference_urls",
  verbose = TRUE
)
```

## Arguments

- data, :

  a [dataset](https://mothur.org/strollur/reference/dataset.md) object

- table, :

  a data.frame containing reference_names, reference_versions
  (optional), reference_usages (optional), reference_notes (optional),
  and reference_urls (optional).

- reference_name, :

  a string containing the name of the column in 'table' that contains
  the reference names. Default column name is 'reference_names'.

- reference_version, :

  a string containing the name of the column in 'table' that contains
  the reference versions. Default column name is 'reference_versions'.

- reference_usage, :

  a string containing the name of the column in 'table' that contains
  the reference usages. Default column name is reference_usages'.

- reference_note, :

  a string containing the name of the column in 'table' that contains
  the reference notes. Default column name is reference_notes'.

- reference_url, :

  a string containing the name of the column in 'table' that contains
  the reference urls. Default column name is reference_urls'.

- verbose, :

  a boolean whether or not you want progress messages. Default = TRUE.

## Value

double containing the number of references added

## Examples

``` r
data <- new_dataset("just for fun", 2)
reference_table <- readr::read_csv(strollur_example("references.csv"),
                             col_names = TRUE, show_col_types = FALSE)
xdev_add_references(data, reference_table)
#> ℹ Added 2 resource references.
#> [1] 2
```
