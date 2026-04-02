# xdev_assign_treatments

Assign samples to treatments in a
[strollur](https://mothur.org/strollur/reference/strollur.md) object

## Usage

``` r
xdev_assign_treatments(
  data,
  table,
  sample = "samples",
  treatment = "treatments",
  verbose = TRUE
)
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.md) object

- table, :

  a data.frame containing sample treatment assignments

- sample, :

  a string containing the name of the column in 'table' that contains
  the samples. Default column name is 'samples'.

- treatment, :

  a string containing the name of the column in 'table' that contains
  the treatments. Default column name is 'treatments'.

- verbose, :

  a boolean indicating whether or not you want progress messages.
  Default = TRUE.

## Value

double containing the number of samples assigned to treatments

## Examples

``` r
data <- new_dataset("my_dataset")
sequence_abundance <- readr::read_tsv(strollur_example(
                                      "mothur2_count_table.tsv.gz"),
                                      show_col_types = FALSE)

xdev_assign_sequence_abundance(data, sequence_abundance, "names")
#> ℹ Assigned 2425 sequence abundances.
#> [1] 2425

sample_assignments <- readr::read_table(file = strollur_example("mouse.time.design"),
                             col_names = TRUE, show_col_types = FALSE)

xdev_assign_treatments(data, sample_assignments)
#> ℹ Assigned 19 samples to treatments.
#> [1] 19
```
