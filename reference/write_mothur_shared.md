# write_mothur_shared

Write mothur formatted [shared
files](https://mothur.org/wiki/shared_file/)

## Usage

``` r
write_mothur_shared(data, file_root = NULL)
```

## Arguments

- data:

  A \`strollur\` object

- file_root:

  a string containing the root name of the output file. Default =
  'dataset_name'. Resulting in output files
  'dataset_name'.bin_type'.shared.

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
write_mothur_shared(miseq, tempfile())
#> [1] "/tmp/Rtmp1sax1K/file1bf9468ca753.otu.shared"      
#> [2] "/tmp/Rtmp1sax1K/file1bf9468ca753.asv.shared"      
#> [3] "/tmp/Rtmp1sax1K/file1bf9468ca753.phylotype.shared"
```
