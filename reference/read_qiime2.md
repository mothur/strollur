# read_qiime2

The read_qiime2 function reads various types of .qza files created by
[qiime2](https://qiime2.org), and creates a \`strollur\` object.

To generate the various input files you can follow [qiime2
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

A \`strollur\` object

## Examples

``` r
# Using the example files from moving-pictures, we add FASTA data, assign
# taxonomy and abundance for features, and add a newick tree and
# metadata.

qza_files <- c(
  strollur_example("rep_seqs.qza"),
  strollur_example("table.qza"),
  strollur_example("taxonomy.qza"),
  strollur_example("rooted-tree.qza")
)

data <- read_qiime2(
  qza = qza_files,
  metadata = strollur_example("sample_metadata.tsv"),
  dataset_name = "qiime2_moving_pictures"
)
#> ℹ Added metadata.
#> ℹ Added 759 sequences.
#> ℹ Assigned 759 sequence abundances.
#> ℹ Assigned 759 asv bins.
#> ℹ Assigned 759 asv bin taxonomies.
data
#> qiime2_moving_pictures:
#> 
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  120    120      0        3     0      1.00
#> 2.5%-tile:       1  120    120      0        3     0   3933.45
#> 25%-tile:        1  120    120      0        4     0  39325.50
#> Median:          1  120    120      0        4     0  78650.00
#> 75%-tile:        1  120    120      0        4     0 117974.50
#> 97.5%-tile:      1  120    120      0        6     0 153366.55
#> Maximum:         1  120    120      0        8     0 157298.00
#> Mean:            1  120    120      0        4     0      0.00
#> 
#> Number of unique seqs: 759 
#> Total number of seqs: 157298 
#> 
#> Total number of samples: 34 
#> Total number of asvs: 759 
#> Your dataset includes metadata 
#> 
```
