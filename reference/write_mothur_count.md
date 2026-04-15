# write_mothur_count

Write a mothur formatted [count
file](https://mothur.org/wiki/count_file/)

## Usage

``` r
write_mothur_count(data, filename = NULL)
```

## Arguments

- data:

  A \`strollur\` object

- filename:

  a string containing the name of the output file. Default =
  'dataset_name'.count_table

## Value

name of count file

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
write_mothur_count(miseq, tempfile())
#> [1] "/tmp/Rtmp1sax1K/file1bf93f84c641"
```
