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

# to create a new dataset that is a copy of miseq

data <- copy_dataset(miseq)
```
