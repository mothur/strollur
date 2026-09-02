# Write mothur formatted [rabund files](https://mothur.org/wiki/rabund_file/)

Write mothur formatted [rabund
files](https://mothur.org/wiki/rabund_file/)

## Usage

``` r
write_mothur_rabund(data, file_root = NULL)
```

## Arguments

- data:

  a [strollur
  object](https://mothur.org/strollur/reference/strollur.html)

- file_root:

  a string containing the root name of the output file. Default =
  'dataset_name'. Resulting in output files
  `{dataset_name}.{bin_type}.rabund`.

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
strollur::write_mothur_rabund(miseq, tempfile())
#> [1] "/tmp/Rtmp5HSFk2/file1d2c3bb7d1d6.otu.rabund"      
#> [2] "/tmp/Rtmp5HSFk2/file1d2c3bb7d1d6.asv.rabund"      
#> [3] "/tmp/Rtmp5HSFk2/file1d2c3bb7d1d6.phylotype.rabund"
```
