# miseq_sop_example

The miseq_sop_example function will create 'dataset' object using the
analysis files from the [MiSeq_SOP](https://mothur.org/wiki/miseq_sop/)
example.

## Usage

``` r
miseq_sop_example()
```

## Value

A 'dataset' object

## Examples

``` r
miseq <- miseq_sop_example()
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
```
