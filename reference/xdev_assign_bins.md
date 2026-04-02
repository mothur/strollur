# xdev_assign_bins

Add bin assignments to a
[strollur](https://mothur.org/strollur/reference/strollur.md) object

## Usage

``` r
xdev_assign_bins(
  data,
  table,
  bin_type = "otu",
  reference = NULL,
  bin_name = "bin_names",
  abundance = "abundances",
  sample = "samples",
  sequence_name = "sequence_names",
  verbose = TRUE
)
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.md) object

- table, :

  a data.frame containing bin_data assignments

- bin_type:

  a string indicating the type of bin assignments. Default "otu".

- reference, :

  a list created by the function \[new_reference\]. Optional.

- bin_name, :

  a string containing the name of the column in 'table' that contains
  the bin names. Default column name is 'bin_names'.

- abundance, :

  a string containing the name of the column in 'table' that contains
  the bin abundances. Default column name is 'abundances'. Note: You
  must provide either abundances or sequence_names in the table.

- sample, :

  a string containing the name of the column in 'table' that contains
  the sample names for datasets where the abundances are broken down by
  sample. Default column name is 'samples'.

- sequence_name, :

  a string containing the name of the column in 'table' that contains
  the sequence names. Default column name is 'sequence_names'. Note: You
  must provide either abundances or sequence_names in the table.

- verbose, :

  a boolean whether or not you want progress messages. Default = TRUE.

## Value

double containing the number of bins assigned

## Examples

``` r
  # To assign sequences to bins:

  data <- new_dataset(dataset_name = "miseq_sop")
  otu_data <- read_mothur_list(list = strollur_example("final.opti_mcc.list.gz"))

  xdev_assign_bins(data = data, table = otu_data, bin_type = "otu")
#> ℹ Assigned 531 otu bins.
#> [1] 531

  # To add abundance only bin assignments:

  data <- new_dataset(dataset_name = "miseq_sop")
  otu_data <- read_mothur_rabund(rabund = strollur_example("final.opti_mcc.rabund"))

  xdev_assign_bins(data = data, table = otu_data, bin_type = "otu")
#> ℹ Assigned 531 otu bins.
#> [1] 531

  # To add abundance bin assignments parsed by sample:

  data <- new_dataset(dataset_name = "miseq_sop")
  otu_data <- readr::read_tsv(strollur_example(
                                "mothur2_bin_assignments_shared.tsv.gz"))
#> Rows: 3431 Columns: 4
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: "\t"
#> chr (3): bin_names, samples, treatments
#> dbl (1): abundances
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

  xdev_assign_bins(data = data, table = otu_data, bin_type = "otu")
#> ℹ Assigned 531 otu bins.
#> [1] 531
```
