# export_dataset

Export all data from a
[dataset](https://mothur.org/strollur/reference/dataset.md) object.

## Usage

``` r
export_dataset(data)
```

## Arguments

- data, :

  a [dataset](https://mothur.org/strollur/reference/dataset.md) object

## Value

Rcpp::List, containing the data in the 'Dataset

## Examples

``` r
dataset <- new_dataset("my_dataset", 2)
export_dataset(dataset)
#> named list()
#> attr(,"strollur_version")
#> [1] "1.0.0"
#> attr(,"dataset_name")
#> [1] "my_dataset"
```
