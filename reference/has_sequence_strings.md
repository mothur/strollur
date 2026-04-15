# has_sequence_strings

Determine if a
[strollur](https://mothur.org/strollur/reference/strollur.html) object
contains sequence nucleotide strings.

## Usage

``` r
has_sequence_strings(data)
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object.

## Value

boolean indicating whether the dataset has sequence nucleotide strings.

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
has_sequence_strings(data)
#> [1] TRUE
```
