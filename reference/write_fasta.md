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
write_fasta(miseq, tempfile())
#> [1] "/tmp/RtmpGZXelw/file1b5d56ccc79f"
```
