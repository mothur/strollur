# export_dataset

Export all data from an instance of the 'Dataset' class.

## Usage

``` r
export_dataset(data, tags = as.character(c()))
```

## Arguments

- data:

  an Rcpp::XPtr\<Dataset\> pointer to an instance of the 'Dataset' c++
  class.

- tags:

  a vector of strings containing the items you wish to export. Options
  are 'sequence_data' and 'bin_data', 'metadata', 'references',
  'sequence_tree', 'sample_tree', and 'reports'. By default, everything
  is exported.

## Value

Rcpp::List, containing the data in the 'Dataset

## Examples

``` r
dataset <- new_dataset("my_dataset", 2)
export_dataset(dataset)
#> named list()
#> attr(,"rdataset_version")
#> [1] "1.0.0"
#> attr(,"dataset_name")
#> [1] "my_dataset"
```
