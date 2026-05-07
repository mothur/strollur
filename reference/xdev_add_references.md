# Add resource references

Add resource references to a
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Usage

``` r
xdev_add_references(
  data,
  table,
  name = "name",
  vendor = "vendor",
  version = "version",
  usage = "usage",
  note = "note",
  method_url = "method_url",
  documentation_url = "documentation_url",
  parameter = "parameter",
  citation = "citation",
  verbose = TRUE
)
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object

- table, :

  a data.frame containing reference_names, reference_versions
  (optional), reference_usages (optional), reference_parameters
  (optional), reference_methods (optional), and reference_urls
  (optional).

- name, :

  a string containing the name of the column in 'table' that contains
  the reference names. Default column name is 'name'.

- vendor, :

  a string containing the name of the column in 'table' that contains
  the reference vendors. Default column name is 'vendor'.

- version, :

  a string containing the name of the column in 'table' that contains
  the reference versions. Default column name is 'version'.

- usage, :

  a string containing the name of the column in 'table' that contains
  the reference usages. Default column name is 'usage'.

- note, :

  a string containing the name of the column in 'table' that contains
  the reference notes. Default column name is 'note'.

- method_url, :

  a string containing the name of the column in 'table' that contains
  the reference methods. Default column name is 'method_url'.

- documentation_url, :

  a string containing the name of the column in 'table' that contains
  the reference documentation urls. Default column name is
  'documentation_url'.

- parameter, :

  a string containing the name of the column in 'table' that contains
  the reference parameters. Default column name is 'parameter'.

- citation, :

  a string containing the name of the column in 'table' that contains
  the reference citations. Default column name is 'citation'.

- verbose, :

  a boolean whether or not you want progress messages. Default = TRUE.

## Value

double containing the number of references added

## Examples

``` r

data <- new_dataset("just for fun")
reference_table <- readr::read_csv(strollur_example("references.csv"),
                             col_names = TRUE, show_col_types = FALSE)
xdev_add_references(data, reference_table)
#> Added 2 resource references.
#> [1] 2
```
