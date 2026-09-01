# is_equal

Determine if two [strollur
object](https://mothur.org/strollur/reference/strollur.html)s are equal.

## Usage

``` r
is_equal(data, data2)
```

## Arguments

- data, :

  a [strollur
  object](https://mothur.org/strollur/reference/strollur.html)

- data2, :

  a [strollur
  object](https://mothur.org/strollur/reference/strollur.html)

## Value

a logical

## Examples

``` r

miseq <- strollur::miseq_sop_example()
#> Added 2425 sequences.
#> Assigned 2425 sequence abundances.
#> Assigned 2425 sequence taxonomies.
#> Assigned 531 otu bins.
#> Assigned 2425 asv bins.
#> Assigned 63 phylotype bins.
#> Assigned 19 samples to treatments.
#> Assigned 171 samples distances.
#> Assigned 531 otu bin taxonomies.
#> Assigned 531 otu bin representative sequences.
#> Added a metadata report.
#> Added 2 resource references.
#> Added a contigs_report report.

data <- strollur::copy_dataset(miseq)

strollur::is_equal(miseq, data)
#> [1] TRUE
```
