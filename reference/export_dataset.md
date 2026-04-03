# export_dataset

Export all data from a
[strollur](https://mothur.org/strollur/reference/strollur.html) object.

## Usage

``` r
export_dataset(data)
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object

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
