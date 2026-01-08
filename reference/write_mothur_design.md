# write_mothur_design

Write a mothur formatted [design
file](https://mothur.org/wiki/design_file/)

## Usage

``` r
write_mothur_design(data, filename = NULL)
```

## Arguments

- data:

  A 'dataset' object

- filename:

  a string containing the name of the output file. Default =
  'dataset_name'.design

## Value

name of design file

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
write_mothur_design(miseq)
#> [1] "miseq_sop.design"
```
