# write_mothur_cons_taxonomy

Write a mothur formatted [cons_taxonomy
file](https://mothur.org/wiki/constaxonomy_file/)

## Usage

``` r
write_mothur_cons_taxonomy(data, file_root = NULL)
```

## Arguments

- data:

  A \`strollur\` object

- file_root:

  a string containing the root name of the output file. Default =
  'dataset_name'. Resulting in output files
  'dataset_name'.bin_type'.cons.taxonomy.

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
write_mothur_cons_taxonomy(miseq, tempfile())
#> [1] "/tmp/RtmpFOmVWd/file1c671b710e7b.otu.cons.taxonomy"      
#> [2] "/tmp/RtmpFOmVWd/file1c671b710e7b.asv.cons.taxonomy"      
#> [3] "/tmp/RtmpFOmVWd/file1c671b710e7b.phylotype.cons.taxonomy"
```
