# Determine if a [strollur object](https://mothur.org/strollur/reference/strollur.html) has sequence taxonomy assignments

Determine if a [strollur
object](https://mothur.org/strollur/reference/strollur.html) has
sequence taxonomy assignments

## Usage

``` r
xdev_has_sequence_taxonomy(data)
```

## Arguments

- data:

  a [strollur
  object](https://mothur.org/strollur/reference/strollur.html)

## Value

boolean

## Examples

``` r

 data <- strollur::miseq_sop_example()
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
 strollur::xdev_has_sequence_taxonomy(data)
#> [1] TRUE
```
