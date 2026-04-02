# xdev_assign_sequence_abundance

Assign sequence abundance and optionally assign sample and treatment
data to a [strollur](https://mothur.org/strollur/reference/strollur.md)
object

## Usage

``` r
xdev_assign_sequence_abundance(
  data,
  table,
  sequence_name = "sequence_names",
  abundance = "abundances",
  sample = "samples",
  treatment = "treatments",
  verbose = TRUE
)
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.md) object

- table, :

  a data.frame containing sequence abundance assignments

- sequence_name, :

  a string containing the name of the column in 'table' that contains
  the sequence names. Default column name is 'sequence_names'.

- abundance, :

  a string containing the name of the column in 'table' that contains
  the sequence abundances. Default column name is 'abundances'.

- sample, :

  a string containing the name of the column in 'table' that contains
  the sequence samples. Default column name is 'samples'. (Optional)

- treatment, :

  a string containing the name of the column in 'table' that contains
  the sequence treatments. Default column name is 'treatments'.

- verbose, :

  a boolean whether or not you want progress messages. Default = TRUE.

## Value

double containing the number of sequences assigned

## Examples

``` r
data <- new_dataset("my_dataset")
sequence_abundance <- readr::read_tsv(strollur_example(
                                      "mothur2_count_table.tsv.gz"),
                                      show_col_types = FALSE)

xdev_assign_sequence_abundance(data = data, table = sequence_abundance,
                               sequence_name = "names")
#> ℹ Assigned 2425 sequence abundances.
#> [1] 2425
```
