# xdev_remove_samples

Designed with package integration in mind, the remove samples function
allows you to remove samples from a [dataset](dataset.md) object

## Usage

``` r
xdev_remove_samples(data, samples)
```

## Arguments

- data, :

  a [dataset](dataset.md) object.

- samples, :

  vector of strings containing the names of the samples to remove.

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

count(data = data, type = "samples")
#> [1] 19

# To remove samples 'F3D0' and 'F3D1'

xdev_remove_samples(data, c("F3D0", "F3D1"))

count(data = data, type = "samples")
#> [1] 17
```
