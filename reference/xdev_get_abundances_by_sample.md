# xdev_get_abundances_by_sample

Get the sequence abundance data in a
[strollur](https://mothur.org/strollur/reference/strollur.html) object
parsed by sample

## Usage

``` r
xdev_get_abundances_by_sample(data, samples = as.character(c()))
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object

- samples:

  a vector of strings containing the names of the samples you would like
  sequence names for. By default all samples are included.

## Value

2D vector of float containing data requested parsed by sample.

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

# To get the sequence names parsed by sample
abunds <- xdev_get_abundances_by_sample(data)
```
