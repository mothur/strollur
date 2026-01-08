# read_qiime2

The read_qiime function reads various types of .qza files created by
[qiime2](https://qiime2.org), and creates a 'dataset' object.

To generate the various input files you can follow [qiime
moving-pictures](https://amplicon-docs.qiime2.org/en/latest/tutorials/moving-pictures.html).

## Usage

``` r
read_qiime2(
  qza,
  metadata = NULL,
  dataset_name = "",
  dir_path = NULL,
  remove_unpacked_artifacts = TRUE
)
```

## Arguments

- qza:

  vector of filenames, .qza files containing your data from qiime2.

- metadata:

  filename, a .tsv file containing metadata

- dataset_name:

  A string containing a name for your dataset.

- dir_path:

  a string containing the name of directory where the artifacts files
  should be unpacked. Default = current working directory.

- remove_unpacked_artifacts:

  boolean, When TRUE, the unpacked artifacts and temporary directories
  will be removed. Default = TRUE.

## Value

A 'dataset' object

## Examples

``` r
# Using the example files from moving-pictures, we add FASTA data, assign
# taxonomy and abundance for features, and add a newick tree and
# metadata.

qza_files <- c(
  rdataset_example("rep_seqs.qza"),
  rdataset_example("table.qza"),
  rdataset_example("taxonomy.qza"),
  rdataset_example("rooted-tree.qza")
)

data <- read_qiime2(
  qza = qza_files,
  metadata = rdataset_example("sample_metadata.tsv"),
  dataset_name = "qiime_moving_pictures"
)
#> ℹ Added metadata.
#> ℹ Added 759 sequences.
#> ℹ Assigned 759 sequence abundances.
#> ℹ Assigned 759 asv bins.
#> ℹ Assigned 759 asv bin taxonomies.
data
#> qiime_moving_pictures:
#> 
#> sequence_summary:
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  120    120      0 3.000000     0      1.00
#> 2.5%-tile:       1  120    120      0 3.000000     0   3933.45
#> 25%-tile:        1  120    120      0 4.000000     0  39325.50
#> Median:          1  120    120      0 4.000000     0  78650.00
#> 75%-tile:        1  120    120      0 4.000000     0 117974.50
#> 97.5%-tile:      1  120    120      0 6.000000     0 153366.55
#> Maximum:         1  120    120      0 8.000000     0 157298.00
#> Mean:            1  120    120      0 4.009059     0      0.00
#> Unique seqs:  759 
#> Total seqs:   157298 
#> 
#> Sample   Total:
#> L1S105   7865 
#> L1S140   7245 
#> L1S208   8270 
#> L1S257   6486 
#> L1S281   6755 
#> L1S57    8756 
#> L1S76    7922 
#> L1S8 7068 
#> L2S155   4112 
#> L2S175   4545 
#> L2S204   3340 
#> L2S222   3485 
#> L2S240   5146 
#> L2S309   1549 
#> L2S357   2526 
#> L2S382   4166 
#> L3S242   917 
#> L3S294   1313 
#> L3S313   1191 
#> L3S341   1109 
#> L3S360   1130 
#> L3S378   1279 
#> L4S112   8575 
#> L4S137   9961 
#> L4S63    10095 
#> L5S104   2253 
#> L5S155   1827 
#> L5S174   1969 
#> L5S203   2132 
#> L5S222   2555 
#> L5S240   1817 
#> L6S20    6892 
#> L6S68    6022 
#> L6S93    7025 
#> 
#> Number of unique seqs: 759 
#> Total number of seqs: 157298 
#> Total number of asvs: 759 
#> 
```
