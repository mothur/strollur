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
#> Added 2425 sequences.
#> Assigned 2425 sequence abundances.
#> Assigned 2425 sequence taxonomies.
#> Assigned 531 otu bins.
#> Assigned 2425 asv bins.
#> Assigned 63 phylotype bins.
#> Assigned 19 samples to treatments.
#> Assigned 531 otu bin taxonomies.
#> Assigned 531 otu bin representative sequences.
#> Added metadata.
#> Added 2 resource references.
#> Added a contigs_report.
 xdev_has_sequence_taxonomy(data)
#> [1] TRUE
```
