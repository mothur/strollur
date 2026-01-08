# has_sample

Determine if a given sample is in a [dataset](dataset.md) object

## Usage

``` r
has_sample(data, sample)
```

## Arguments

- data, :

  a [dataset](dataset.md) object.

- sample:

  a string containing the name of a sample.

## Value

boolean indicating whether the dataset has a given sample

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
has_sample(data, "F3D0")
#> [1] TRUE
has_sample(data, "not a valid sample")
#> [1] FALSE
```
