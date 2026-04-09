# write_taxonomy

Write a 2 column [taxonomy file](https://mothur.org/wiki/taxonomy_file/)

## Usage

``` r
write_taxonomy(data, filename = NULL)
```

## Arguments

- data:

  A \`strollur\` object

- filename:

  a string containing the name of the output file. Default =
  'dataset_name'.taxonomy

## Value

name of taxonomy file

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
write_taxonomy(miseq, tempfile())
#> [1] "/tmp/RtmpOaitIa/file1bcc244eb531"
```
