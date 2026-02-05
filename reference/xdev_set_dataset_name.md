# xdev_set_dataset_name

Designed with package integration in mind, set the name of a
[dataset](https://mothur.org/strollur/reference/dataset.md) object.

## Usage

``` r
xdev_set_dataset_name(data, dataset_name)
```

## Arguments

- data, :

  a [dataset](https://mothur.org/strollur/reference/dataset.md) object

- dataset_name, :

  a string containing the desired name

## Examples

``` r
data <- new_dataset(dataset_name = "my_dataset")
xdev_set_dataset_name(data = data, dataset_name = "new_dataset_name")
```
