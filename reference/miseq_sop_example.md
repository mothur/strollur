# Create a [strollur object](https://mothur.org/strollur/reference/strollur.html) using the analysis files from the [MiSeq_SOP](https://mothur.org/wiki/miseq_sop/) example.

The miseq_sop_example function will create
[strollur](https://mothur.org/strollur/reference/strollur.html) object
using the analysis files from the
[MiSeq_SOP](https://mothur.org/wiki/miseq_sop/) example.

## Usage

``` r
miseq_sop_example()
```

## Value

A
[`strollur::strollur`](https://mothur.org/strollur/reference/strollur.md)
object

## Examples

``` r

miseq <- strollur::miseq_sop_example()
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
```
