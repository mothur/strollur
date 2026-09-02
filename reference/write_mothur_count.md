# Write a mothur formatted [count file](https://mothur.org/wiki/count_file/)

Write a mothur formatted [count
file](https://mothur.org/wiki/count_file/)

## Usage

``` r
write_mothur_count(data, filename = NULL)
```

## Arguments

- data:

  a [strollur
  object](https://mothur.org/strollur/reference/strollur.html)

- filename:

  a string containing the name of the output file. Default =
  'dataset_name'.count_table

## Value

name of count file

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
strollur::write_mothur_count(miseq, tempfile())
#> [1] "/tmp/Rtmp5HSFk2/file1d2c6d8f1aa0"
```
