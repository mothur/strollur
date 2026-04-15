# write_fasta

Write a [FASTA](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/)
formatted sequence file

## Usage

``` r
write_fasta(data, filename = NULL, degap = FALSE)
```

## Arguments

- data:

  A \`strollur\` object

- filename:

  a string containing the name of the output file. Default =
  'dataset_name'.fasta

- degap:

  a logical. Default = FALSE. When degap = \`TRUE\`, all gap characters
  will be removed from the sequences.

## Value

name of FASTA file

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
write_fasta(miseq, tempfile())
#> [1] "/tmp/Rtmp1sax1K/file1bf938f64d2d"
```
