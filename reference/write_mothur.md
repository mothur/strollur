# write_mothur

The write_mothur function will write various [file
types](https://mothur.org/wiki/tags/#file_types) for use with mothur.

## Usage

``` r
write_mothur(data, dir_path = NULL, compress = TRUE, tags = NULL)
```

## Arguments

- data:

  A 'dataset' object

- dir_path:

  a string containing the name of directory where the files should be
  written. Default = current working directory.

- compress:

  boolean, Default = TRUE.

- tags:

  a vector of strings containing the items you wish to write Options are
  'sequence_data', 'bin_data', 'metadata', 'references',
  'sequence_tree', 'sample_tree' and 'reports'. By default, everything
  is written to files.

## Value

a vector of file names

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
files <- write_mothur(miseq)
```
