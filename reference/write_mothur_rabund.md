# write_mothur_rabund

Write mothur formatted [rabund
files](https://mothur.org/wiki/rabund_file/)

## Usage

``` r
write_mothur_rabund(data, file_root = NULL)
```

## Arguments

- data:

  A \`strollur\` object

- file_root:

  a string containing the root name of the output file. Default =
  'dataset_name'. Resulting in output files
  'dataset_name'.bin_type'.rabund.

## Value

vector containing the names of the files created

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
write_mothur_rabund(miseq, tempfile())
#> [1] "/tmp/Rtmp5CtWbg/file1b5428d4aa8c.otu.rabund"      
#> [2] "/tmp/Rtmp5CtWbg/file1b5428d4aa8c.asv.rabund"      
#> [3] "/tmp/Rtmp5CtWbg/file1b5428d4aa8c.phylotype.rabund"
```
