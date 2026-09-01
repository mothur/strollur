# Create a new [strollur object](https://mothur.org/strollur/reference/strollur.html)

Create a new [strollur
object](https://mothur.org/strollur/reference/strollur.html)

## Usage

``` r
new_dataset(dataset_name = "")
```

## Arguments

- dataset_name:

  string, a string containing the dataset name. Default = ""

## Value

a [strollur object](https://mothur.org/strollur/reference/strollur.html)

## See also

The 'new' method in the [strollur
object](https://mothur.org/strollur/reference/strollur.html)

## Examples

``` r

data <- strollur::new_dataset()

# to create a new dataset named "soil", run the following:

data <- strollur::new_dataset(dataset_name = "soil")
```
