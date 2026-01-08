# write_mothur_list

Write mothur formatted [list files](https://mothur.org/wiki/list_file/)

## Usage

``` r
write_mothur_list(data, file_root = NULL)
```

## Arguments

- data:

  A 'dataset' object

- file_root:

  a string containing the root name of the output file. Default =
  'dataset_name'. Resulting in output files
  'dataset_name'.bin_type'.list.

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
write_mothur_list(miseq)
#> [1] "miseq_sop.otu.list"       "miseq_sop.asv.list"      
#> [3] "miseq_sop.phylotype.list"
```
