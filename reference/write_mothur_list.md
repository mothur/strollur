# write_mothur_list

Write mothur formatted [list files](https://mothur.org/wiki/list_file/)

## Usage

``` r
write_mothur_list(data, file_root = NULL)
```

## Arguments

- data:

  A \`strollur\` object

- file_root:

  a string containing the root name of the output file. Default =
  'dataset_name'. Resulting in output files
  'dataset_name'.bin_type'.list.

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
#> Added a metadata report.
#> Added 2 resource references.
#> Added a contigs_report report.
write_mothur_list(miseq, tempfile())
#> [1] "/tmp/RtmpHbXlcd/file1b0a2e60363b.otu.list"      
#> [2] "/tmp/RtmpHbXlcd/file1b0a2e60363b.asv.list"      
#> [3] "/tmp/RtmpHbXlcd/file1b0a2e60363b.phylotype.list"
```
