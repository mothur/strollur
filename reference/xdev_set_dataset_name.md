# Set the name of a [strollur object](https://mothur.org/strollur/reference/strollur.html)

Designed with package integration in mind, set the name of a [strollur
object](https://mothur.org/strollur/reference/strollur.html).

## Usage

``` r
xdev_set_dataset_name(data, dataset_name)
```

## Arguments

- data:

  a [strollur
  object](https://mothur.org/strollur/reference/strollur.html)

- dataset_name:

  a string containing the desired name

## Value

No return value, called for side effects.

## Examples

``` r

data <- strollur::new_dataset(dataset_name = "my_dataset")
strollur::xdev_set_dataset_name(data,
                                dataset_name = "new_dataset_name")
```
