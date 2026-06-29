# is_equal

Determine if two
[strollur](https://mothur.org/strollur/reference/strollur.html) objects
are equal.

## Usage

``` r
is_equal(data, data2)
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object

- data2, :

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object

## Value

a logical

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
#> Added a metadata report.
#> Added 2 resource references.
#> Added a contigs_report report.

data <- copy_dataset(miseq)

is_equal(miseq, data)
#> [1] TRUE
```
