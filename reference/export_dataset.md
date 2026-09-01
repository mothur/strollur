# Create a human readable table containing all data from a [strollur](https://mothur.org/strollur/reference/strollur.html) object.

Export all data from a
[strollur](https://mothur.org/strollur/reference/strollur.html) object.

## Usage

``` r
export_dataset(data)
```

## Arguments

- data:

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object

## Value

list, containing the data in the strollur::strollur object

## Examples

``` r

dataset <- strollur::new_dataset("my_dataset")
strollur::export_dataset(dataset)
#> named list()
#> attr(,"strollur_version")
#> [1] "0.1.2"
#> attr(,"dataset_name")
#> [1] "my_dataset"
```
