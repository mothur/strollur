# new_dataset

Create a new
[strollur](https://mothur.org/strollur/reference/strollur.md) object

## Usage

``` r
new_dataset(dataset_name = "", processors = NULL)
```

## Arguments

- dataset_name:

  string, a string containing the dataset name. Default = ""

- processors:

  integer, number of cores to use during summary functions. Default =
  all available

## Value

a [strollur](https://mothur.org/strollur/reference/strollur.md) object

## See also

The 'new' method in the
[strollur](https://mothur.org/strollur/reference/strollur.md) class

## Examples

``` r
data <- new_dataset()

# to create a new dataset named "soil" and allow for all available
# processors during summary functions, run the following:

data <- new_dataset(dataset_name = "soil")

# to create a new dataset named "soil" and allow for 2
# processors during summary functions, run the following:

data <- new_dataset(dataset_name = "soil", processors = 2)
```
