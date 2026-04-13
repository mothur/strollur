# assigned_message

Report dataset assignments

## Usage

``` r
assigned_message(num, tag = "sequences")
```

## Arguments

- num:

  integer containing the number of items assigned

- tag:

  string containing item assigned Default = 'sequences'.

## Value

No return value, called for side effects.

## Examples

``` r
assigned_message(10)
#> ℹ Assigned 10sequences
assigned_message(10, "sequence abundances")
#> ℹ Assigned 10sequence abundances
```
