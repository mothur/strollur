# is_aligned

Determine if a [dataset](dataset.md) object contains aligned sequences.

## Usage

``` r
is_aligned(data)
```

## Arguments

- data, :

  a [dataset](dataset.md) object

## Value

Boolean

## Examples

``` r
dataset <- miseq_sop_example()
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
is_aligned(dataset)
#> [1] TRUE
```
