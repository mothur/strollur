# xdev_assign_bin_taxonomy

Assign bin classifications to a
[dataset](https://mothur.org/strollur/reference/dataset.md) object

Note, if you assign sequence taxonomies and assign bins, 'Dataset' will
find the concensus taxonomy for each bin for you.

## Usage

``` r
xdev_assign_bin_taxonomy(
  data,
  table,
  bin_type = "otu",
  reference = NULL,
  bin_name = "bin_names",
  taxonomy = "taxonomies",
  verbose = TRUE
)
```

## Arguments

- data, :

  a [dataset](https://mothur.org/strollur/reference/dataset.md) object

- table, :

  a data.frame containing bin taxonomy assignments

- bin_type:

  a string indicating the type of bin assignments. Default "otu".

- reference, :

  a list created by the function \[new_reference\]. Optional.

- bin_name, :

  a string containing the name of the column in 'table' that contains
  the bin names. Default column name is 'bin_names'.

- taxonomy, :

  a string containing the name of the column in 'table' that contains
  the bin taxonomies. Default column name is 'taxonomies'.

- verbose, :

  a boolean whether or not you want progress messages. Default = TRUE.

## Value

double containing the number of bins assigned

## Examples

``` r
otu_data <- read_mothur_cons_taxonomy(strollur_example(
                        "final.cons.taxonomy"))

data <- new_dataset(dataset_name = "my_dataset")

# assign otu abundances
xdev_assign_bins(data = data, table = otu_data, bin_type = "otu")
#> ℹ Assigned 531 otu bins.
#> [1] 531

# assign otu classifications
xdev_assign_bin_taxonomy(data = data, table = otu_data,
                         bin_type = "otu")
#> ℹ Assigned 531 otu bin taxonomies.
#> [1] 531
```
