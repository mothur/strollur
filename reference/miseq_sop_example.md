# Example strollur object

The miseq_sop_example function will create 'strollur' object using the
analysis files from the [MiSeq_SOP](https://mothur.org/wiki/miseq_sop/)
example.

## Usage

``` r
miseq_sop_example()
```

## Value

A 'strollur' object

## Examples

``` r
miseq <- miseq_sop_example()
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
```
