# xdev_remove_samples

Designed with package integration in mind, the remove samples function
allows you to remove samples from a
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Usage

``` r
xdev_remove_samples(data, samples, reason = "remove_samples")
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object.

- samples, :

  vector of strings containing the names of the samples to remove.

- reason, :

  string containing the reason for removal. Default = "remove_samples".

## Value

No return value, called for side effects.

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

count(data = data, type = "samples")
#> [1] 19

# To remove samples 'F3D0' and 'F3D1'

xdev_remove_samples(data, c("F3D0", "F3D1"))

count(data = data, type = "samples")
#> [1] 17
```
