# xdev_set_num_processors

Designed with package integration in mind, set the number of processors
used to summarize a
[strollur](https://mothur.org/strollur/reference/strollur.md) object

## Usage

``` r
xdev_set_num_processors(data, processors)
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.md) object

- processors, :

  a integer containing the desired number of processors

## Examples

``` r
data <- new_dataset(dataset_name = "my_dataset")
xdev_set_num_processors(data = data, processors = 1)
```
