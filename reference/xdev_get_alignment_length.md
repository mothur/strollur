# Get the alignment length of sequences in your [strollur object](https://mothur.org/strollur/reference/strollur.html)

Get the alignment length of sequences in your [strollur
object](https://mothur.org/strollur/reference/strollur.html)

## Usage

``` r
xdev_get_alignment_length(data)
```

## Arguments

- data:

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object

## Value

Integer containing the length of the alignment or -1 if unaligned.

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
 strollur::xdev_get_alignment_length(data)
#> [1] 375
```
