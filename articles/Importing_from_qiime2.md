# Importing from qiime2

*strollur* includes the function
[`read_qiime2()`](https://mothur.org/strollur/reference/read_qiime2.md)
as well as several functions to read [qiime2](https://qiime2.org) output
files individually. To create a dataset from the outputs of the [qiime2
moving-pictures
example](https://amplicon-docs.qiime2.org/en/latest/tutorials/moving-pictures.html),
run the following:

``` r
qza_files <- c(
  strollur_example("rep_seqs.qza"),
  strollur_example("table.qza"),
  strollur_example("taxonomy.qza"),
  strollur_example("rooted-tree.qza")
)

data <- read_qiime2(
  qza = qza_files,
  metadata = strollur_example("sample_metadata.tsv"),
  dataset_name = "qiime_moving_pictures"
)
#> ℹ Added metadata.
#> ℹ Added 759 sequences.
#> ℹ Assigned 759 sequence abundances.
#> ℹ Assigned 759 asv bins.
#> ℹ Assigned 759 asv bin taxonomies.
```

To view a summary of data:

``` r
data
#> qiime_moving_pictures:
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
#> Total number of asv bin classifications: 759 
#> Your dataset includes metadata
```

## Reading Individual Files

- [`unpack_qiime2_artifact()`](https://mothur.org/strollur/reference/unpack_qiime2_artifact.md)
  unpacks qiime2 *qza* files and returns an artifact
- [`read_qiime2_feature_table()`](https://mothur.org/strollur/reference/read_qiime2_feature_table.md)
  reads a qiime2 *qza* file containing bin data
- [`read_qiime2_taxonomy()`](https://mothur.org/strollur/reference/read_qiime2_taxonomy.md)
  reads a qiime2 *qza* file containing containing taxonomy data
- [`read_qiime2_metadata()`](https://mothur.org/strollur/reference/read_qiime2_metadata.md)
  read qiime2 *tsv* table containing metadata

To create a dataset and read the individual files, you can use the
functions below. First let’s create a dataset named my_data.

``` r
my_data <- new_dataset(dataset_name = "my_data")
```

To decompress the individual qza files, you can use the
[`unpack_qiime2_artifact()`](https://mothur.org/strollur/reference/unpack_qiime2_artifact.md).
unpack_qiime2_artifact will extract a folder containing the data,
versioning and provenance information. Let’s take a look at the
rep_seqs.qza file first.

``` r
fasta_artifact <- unpack_qiime2_artifact(qza = strollur_example("rep_seqs.qza"))
```

fasta_artifact is a list. To add the
[FASTA](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/) data to your
dataset you can use the
[`read_fasta()`](https://mothur.org/strollur/reference/read_fasta.md)
function:

``` r
fasta_file <- file.path(
  getwd(), "rep_seqs", fasta_artifact$uuid, "data",
  "dna-sequences.fasta"
)
fasta_data <- read_fasta(fasta = fasta_file)
```

fasta_data is a data.frame containing sequence names, sequence
nucleotide strings, and comments if provided. You can add the FASTA
sequences to your dataset using the
[`add()`](https://mothur.org/strollur/reference/add.md) function:

``` r
add(
  my_data,
  table = fasta_data,
  type = "sequences"
)
#> ℹ Added 759 sequences.
#> [1] 759
```

You can remove the fasta_artifact using the unlink function:

``` r
unlink(file.path(getwd(), "rep_seqs"), recursive = TRUE)
```

Now that we have added the FASTA sequences, let’s read the abundance
data found in the table.qza file. Since we know that the sequences are
in fact features, we will assign the abundances to the sequences and
then we can assign the sequences to bins to create [Amplicon Sequence
Variant](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) *asv*
clusters. Let’s use the
[`read_qiime2_feature_table()`](https://mothur.org/strollur/reference/read_qiime2_feature_table.md)
function to extract the abundance table.

``` r
abundance <- read_qiime2_feature_table(strollur_example("table.qza"))

assign(
  my_data,
  table = abundance$data,
  type = "sequence_abundance",
  table_names = list(sequence_name = "bin_names")
)
#> ℹ Assigned 759 sequence abundances.
#> [1] 759

assign(
  my_data,
  table = abundance$data,
  type = "bins",
  bin_type = "asv",
  table_names = list(sequence_name = "bin_names")
)
#> ℹ Assigned 759 asv bins.
#> [1] 759

my_data
#> my_data:
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
```

Now that we have added the abundance data and assigned the sequences to
bins, let’s take a look at the taxonomic classifications provided in the
taxonomy.qza file. We will use the
[`read_qiime2_taxonomy()`](https://mothur.org/strollur/reference/read_qiime2_taxonomy.md)
function to extract the feature classification table.

``` r
taxonomy <- read_qiime2_taxonomy(strollur_example("taxonomy.qza"))

assign(
  my_data,
  table = taxonomy$data,
  type = "sequence_taxonomy",
  table_names = list(sequence_name = "bin_names")
)
#> ℹ Assigned 759 sequence taxonomies.
#> [1] 759
```

Now, lets add a tree that shows the relationships between the sequences
(features).

``` r
tree_artifact <- unpack_qiime2_artifact(
  qza =
    strollur_example("rooted-tree.qza")
)
```

We can read the tree using the
[`ape::read.tree()`](https://rdrr.io/pkg/ape/man/read.tree.html)
function.

``` r
tree_file <- file.path(
  getwd(), "rooted-tree", tree_artifact$uuid, "data",
  "tree.nwk"
)
sequence_tree <- ape::read.tree(tree_file)
```

You can add the tree to your dataset as follows:

``` r
my_data$add_sequence_tree(tree = sequence_tree)
```

You can remove the tree_artifact using the following:

``` r
unlink(file.path(getwd(), "rooted-tree"), recursive = TRUE)
```

Lastly, let’s read the metadata provided by qiime2 using the
[`read_qiime2_metadata()`](https://mothur.org/strollur/reference/read_qiime2_metadata.md)
function.

``` r
metadata <- read_qiime2_metadata(
  metadata =
    strollur_example("sample_metadata.tsv")
)
```

To add the metadata to the my_data dataset, run the following:

``` r
add(my_data, table = metadata, type = "metadata")
#> ℹ Added metadata.
#> [1] 1
```

To view a summary of your imported dataset, ‘my_data’, run the
following:

``` r
my_data
#> my_data:
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
#> Total number of asv bin classifications: 759 
#> Total number of sequence classifications: 759 
#> Your dataset includes metadata
```
