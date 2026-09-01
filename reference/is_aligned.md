# Determine if a [strollur object](https://mothur.org/strollur/reference/strollur.html) contains aligned sequences.

Determine if a [strollur
object](https://mothur.org/strollur/reference/strollur.html) contains
aligned sequences.

## Usage

``` r
is_aligned(data)
```

## Arguments

- data, :

  a [strollur
  object](https://mothur.org/strollur/reference/strollur.html)

## Value

Boolean

## Examples

``` r

dataset <- miseq_sop_example()
#> Added 2425 sequences.
#> Assigned 2425 sequence abundances.
#> Assigned 2425 sequence taxonomies.
#> Assigned 531 otu bins.
#> Assigned 2425 asv bins.
#> Assigned 63 phylotype bins.
#> Assigned 19 samples to treatments.
#> Assigned 171 samples distances.
#> Assigned 531 otu bin taxonomies.
#> Assigned 531 otu bin representative sequences.
#> Added a metadata report.
#> Added 2 resource references.
#> Added a contigs_report report.
is_aligned(dataset)
#> [1] TRUE
```
