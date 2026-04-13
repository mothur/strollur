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
write_mothur_shared(miseq, tempfile())
#> [1] "/tmp/RtmpGZXelw/file1b5d424acc5d.otu.shared"      
#> [2] "/tmp/RtmpGZXelw/file1b5d424acc5d.asv.shared"      
#> [3] "/tmp/RtmpGZXelw/file1b5d424acc5d.phylotype.shared"
```
