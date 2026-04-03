# copy_dataset

Create a new
[strollur](https://mothur.org/strollur/reference/strollur.html) object
from an existing dataset.

## Usage

``` r
copy_dataset(data)
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object

## Value

a [strollur](https://mothur.org/strollur/reference/strollur.html) object

## See also

The 'new' method in the
[strollur](https://mothur.org/strollur/reference/strollur.html) class

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

# to create a new dataset that is a copy of miseq

data <- copy_dataset(miseq)
```
