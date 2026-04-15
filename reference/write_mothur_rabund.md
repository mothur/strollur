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
write_mothur_rabund(miseq, tempfile())
#> [1] "/tmp/Rtmp1sax1K/file1bf91f46404d.otu.rabund"      
#> [2] "/tmp/Rtmp1sax1K/file1bf91f46404d.asv.rabund"      
#> [3] "/tmp/Rtmp1sax1K/file1bf91f46404d.phylotype.rabund"
```
