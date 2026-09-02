# Write mothur formatted [list files](https://mothur.org/wiki/list_file/)

Write mothur formatted [list files](https://mothur.org/wiki/list_file/)

## Usage

``` r
write_mothur_list(data, file_root = NULL)
```

## Arguments

- data:

  a [strollur
  object](https://mothur.org/strollur/reference/strollur.html)

- file_root:

  a string containing the root name of the output file. Default =
  'dataset_name'. Resulting in output files
  `{dataset_name}.{bin_type}.list`.

## Value

vector containing the names of the files created

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
strollur::write_mothur_list(miseq, tempfile())
#> [1] "/tmp/Rtmp5HSFk2/file1d2c36597cc7.otu.list"      
#> [2] "/tmp/Rtmp5HSFk2/file1d2c36597cc7.asv.list"      
#> [3] "/tmp/Rtmp5HSFk2/file1d2c36597cc7.phylotype.list"
```
