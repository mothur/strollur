# xdev_set_dataset_name

Designed with package integration in mind, set the name of a
[strollur](https://mothur.org/strollur/reference/strollur.html) object.

## Usage

``` r
xdev_set_dataset_name(data, dataset_name)
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object

- dataset_name, :

  a string containing the desired name

## Value

No return value, called for side effects.

## Examples

``` r
data <- new_dataset(dataset_name = "my_dataset")
xdev_set_dataset_name(data = data, dataset_name = "new_dataset_name")
```
