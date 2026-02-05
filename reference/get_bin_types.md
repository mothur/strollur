# get_bin_types

Get bin table types of a
[dataset](https://mothur.org/strollur/reference/dataset.md) object

## Usage

``` r
get_bin_types(data)
```

## Arguments

- data, :

  a [dataset](https://mothur.org/strollur/reference/dataset.md) object

## Value

vector of strings

## Examples

``` r
data <- miseq_sop_example()
#> ℹ Added 2425 sequences.
#> ℹ Assigned 2425 sequence abundances.
#> ℹ Assigned 2425 sequence taxonomies.
#> ℹ Assigned 531 otu bins.
#> ℹ Assigned 2425 asv bins.
#> ℹ Assigned 63 phylotype bins.
#> ℹ Assigned 19 samples to treatments.
#> ℹ Assigned 531 otu bin taxonomies.
#> ℹ Assigned 531 otu bin representative sequences.
#> ℹ Added metadata.
#> ℹ Added 2 resource references.
#> ℹ Added a contigs_report.
get_bin_types(data)
#> [1] "otu"       "asv"       "phylotype"
```
