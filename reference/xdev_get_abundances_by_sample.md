# xdev_get_abundances_by_sample

Get the sequence abundance data in a
[strollur](https://mothur.org/strollur/reference/strollur.md) object
parsed by sample

## Usage

``` r
xdev_get_abundances_by_sample(data, samples = as.character(c()))
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.md) object

- samples:

  a vector of strings containing the names of the samples you would like
  sequence names for. By default all samples are included.

## Value

2D vector of float containing data requested parsed by sample.

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

# To get the sequence names parsed by sample
abunds <- xdev_get_abundances_by_sample(data)
```
