# xdev_has_sequence_taxonomy

Determine if a
[strollur](https://mothur.org/strollur/reference/strollur.html) object
has sequence taxonomy assignments

## Usage

``` r
xdev_has_sequence_taxonomy(data)
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object

## Value

boolean

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
 xdev_has_sequence_taxonomy(data)
#> [1] TRUE
```
