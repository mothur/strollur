# added_message

Report dataset additions

## Usage

``` r
added_message(num, tag = "sequences")
```

## Arguments

- num:

  integer containing the number of items added.

- tag:

  string containing item added. Default = 'sequences'.

## Value

No return value, called for side effects.

## Examples

``` r
added_message(10)
#> ℹ Added 10 sequences.
added_message(10, "sequence abundances")
#> ℹ Added 10 sequence abundances.
```
