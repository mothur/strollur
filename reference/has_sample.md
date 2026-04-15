# has_sample

Determine if a given sample is in a
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Usage

``` r
has_sample(data, sample)
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object.

- sample:

  a string containing the name of a sample.

## Value

boolean indicating whether the dataset has a given sample

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
has_sample(data, "F3D0")
#> [1] TRUE
has_sample(data, "not a valid sample")
#> [1] FALSE
```
