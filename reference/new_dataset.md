# new_dataset

Create a new
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Usage

``` r
new_dataset(dataset_name = "")
```

## Arguments

- dataset_name:

  string, a string containing the dataset name. Default = ""

## Value

a [strollur](https://mothur.org/strollur/reference/strollur.html) object

## See also

The 'new' method in the
[strollur](https://mothur.org/strollur/reference/strollur.html) class

## Examples

``` r

data <- new_dataset()

# to create a new dataset named "soil", run the following:

data <- new_dataset(dataset_name = "soil")
```
