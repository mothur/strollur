# dataset

'dataset' is an R6 class that stores nucleotide sequences, abundance,
sample and treatment assignments, taxonomic classifications, asv / otu
clusters and various reports. It is designed to facilitate data analysis
across multiple R packages.

## Author

Sarah Westcott, <swestcot@umich.edu>

## Public fields

- `data`:

  Rcpp::XPtr\<Dataset\> pointer to 'Dataset' c++ class. This allows
  package developers an easy access point to the underlying C++ code
  with additional functionality.

- `raw`:

  Rcpp::RawVector containing the serialized data of the 'Dataset' c++
  class. This allows the load and save functions to work with the class.

- `sequence_tree`:

  a tree that relates sequences to eachother

- `sample_tree`:

  a tree that relates samples to eachother

## Methods

### Public methods

- [`dataset$new()`](#method-dataset-new)

- [`dataset$print()`](#method-dataset-print)

- [`dataset$add_sample_tree()`](#method-dataset-add_sample_tree)

- [`dataset$add_sequence_tree()`](#method-dataset-add_sequence_tree)

- [`dataset$clear()`](#method-dataset-clear)

- [`dataset$get_bin_types()`](#method-dataset-get_bin_types)

- [`dataset$get_metadata()`](#method-dataset-get_metadata)

- [`dataset$get_sample_tree()`](#method-dataset-get_sample_tree)

- [`dataset$get_sequence_report()`](#method-dataset-get_sequence_report)

- [`dataset$get_summary()`](#method-dataset-get_summary)

- [`dataset$get_sequence_tree()`](#method-dataset-get_sequence_tree)

- [`dataset$clone()`](#method-dataset-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new dataset

#### Usage

    dataset$new(
      name = "",
      processors = parallelly::availableCores(),
      dataset = NULL
    )

#### Arguments

- `name`:

  String, name of dataset (optional)

- `processors`:

  Integer, number of cores to use. Default = all available

- `dataset`:

  a \`dataset\` object.

#### Returns

A new \`dataset\` object.

#### Examples

    # to create an empty dataset, run the following:

    data <- new_dataset("soil")

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Get summary of sequence data

#### Usage

    dataset$print()

------------------------------------------------------------------------

### Method `add_sample_tree()`

Add phylo tree relating the samples in your dataset

#### Usage

    dataset$add_sample_tree(tree)

#### Arguments

- `tree`:

  a phylo tree object created by ape::read.tree.

#### Examples

     data <- dataset$new("my_dataset")

     df <- read_mothur_shared(strollur_example("final.opti_mcc.shared"))
     assign(data = data, table = df, type = "bins", bin_type = "otu")

     tree <- ape::read.tree(strollur_example(
     "final.opti_mcc.jclass.ave.tre"))

     data$add_sample_tree(tree)

------------------------------------------------------------------------

### Method `add_sequence_tree()`

Add phylo tree relating the sequences in your dataset

#### Usage

    dataset$add_sequence_tree(tree)

#### Arguments

- `tree`:

  a phylo tree object created by ape::read.tree.

#### Examples

     data <- dataset$new("my_dataset")
     tree <- ape::read.tree(strollur_example("final.phylip.tre.gz"))
     data$add_sequence_tree(tree)

------------------------------------------------------------------------

### Method [`clear()`](https://mothur.org/strollur/reference/clear.md)

Clear data from datasest

#### Usage

    dataset$clear()

------------------------------------------------------------------------

### Method [`get_bin_types()`](https://mothur.org/strollur/reference/get_bin_types.md)

Get bin table types

#### Usage

    dataset$get_bin_types()

#### Returns

vector of strings

#### Examples

    data <- miseq_sop_example()
    data$get_bin_types()

------------------------------------------------------------------------

### Method `get_metadata()`

Get data.frame containing metadata for the dataset

#### Usage

    dataset$get_metadata()

#### Returns

data.frame()

#### Examples

      data <- dataset$new("my_dataset")

      metadata <- readr::read_tsv(strollur_example("sample-metadata.tsv"),
       col_names = TRUE, show_col_types = FALSE)

      add(data = data, table = metadata, type = "metadata")

      data$get_metadata()

------------------------------------------------------------------------

### Method `get_sample_tree()`

Get phylo tree relating the samples in your dataset.

#### Usage

    dataset$get_sample_tree()

#### Examples

     tree <- ape::read.tree(strollur_example(
      "final.opti_mcc.jclass.ave.tre"))

     df <- read_mothur_shared(strollur_example("final.opti_mcc.shared"))

     data <- dataset$new("my_dataset")

     # assign abundance 'otu' bins
     assign(data = data, table = df, type = "bins", bin_type = "otu")

     data$add_sample_tree(tree)
     data$get_sample_tree()

------------------------------------------------------------------------

### Method `get_sequence_report()`

Get data.frame sequence report data. Sequence report data includes:
start positions, end positions, number of bases, number of ambiguous
bases, length of longest homopolymer, and the number of N's.

#### Usage

    dataset$get_sequence_report()

#### Returns

data.frame

------------------------------------------------------------------------

### Method `get_summary()`

Get summary of the sequence reports

#### Usage

    dataset$get_summary(silent = FALSE)

#### Arguments

- `silent`:

  Default = FALSE, meaning print summaries

#### Returns

list of data.frames

------------------------------------------------------------------------

### Method `get_sequence_tree()`

Get phylo tree relating the sequences in your dataset.

#### Usage

    dataset$get_sequence_tree()

#### Examples

     data <- dataset$new("my_dataset")
     tree <- ape::read.tree(strollur_example("final.phylip.tre.gz"))
     data$add_sequence_tree(tree)
     data$get_sequence_tree()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    dataset$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
## ------------------------------------------------
## Method `dataset$new`
## ------------------------------------------------


# to create an empty dataset, run the following:

data <- new_dataset("soil")


## ------------------------------------------------
## Method `dataset$add_sample_tree`
## ------------------------------------------------


 data <- dataset$new("my_dataset")

 df <- read_mothur_shared(strollur_example("final.opti_mcc.shared"))
 assign(data = data, table = df, type = "bins", bin_type = "otu")
#> ℹ Assigned 531 otu bins.
#> [1] 531

 tree <- ape::read.tree(strollur_example(
 "final.opti_mcc.jclass.ave.tre"))

 data$add_sample_tree(tree)


## ------------------------------------------------
## Method `dataset$add_sequence_tree`
## ------------------------------------------------


 data <- dataset$new("my_dataset")
 tree <- ape::read.tree(strollur_example("final.phylip.tre.gz"))
 data$add_sequence_tree(tree)
#> ℹ Added 2425 sequences.


## ------------------------------------------------
## Method `dataset$get_bin_types`
## ------------------------------------------------


data <- miseq_sop_example()
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
data$get_bin_types()
#> [1] "otu"       "asv"       "phylotype"


## ------------------------------------------------
## Method `dataset$get_metadata`
## ------------------------------------------------

  data <- dataset$new("my_dataset")

  metadata <- readr::read_tsv(strollur_example("sample-metadata.tsv"),
   col_names = TRUE, show_col_types = FALSE)

  add(data = data, table = metadata, type = "metadata")
#> ℹ Added metadata.
#> [1] 1

  data$get_metadata()
#>    sample-id barcode-sequence  body-site year month day   subject
#> 1       L1S8     AGCTGACTAGTC        gut 2008    10  28 subject-1
#> 2      L1S57     ACACACTATGGC        gut 2009     1  20 subject-1
#> 3      L1S76     ACTACGTGTGGT        gut 2009     2  17 subject-1
#> 4     L1S105     AGTGCGATGCGT        gut 2009     3  17 subject-1
#> 5     L2S155     ACGATGCGACCA  left palm 2009     1  20 subject-1
#> 6     L2S175     AGCTATCCACGA  left palm 2009     2  17 subject-1
#> 7     L2S204     ATGCAGCTCAGT  left palm 2009     3  17 subject-1
#> 8     L2S222     CACGTGACATGT  left palm 2009     4  14 subject-1
#> 9     L3S242     ACAGTTGCGCGA right palm 2008    10  28 subject-1
#> 10    L3S294     CACGACAGGCTA right palm 2009     1  20 subject-1
#> 11    L3S313     AGTGTCACGGTG right palm 2009     2  17 subject-1
#> 12    L3S341     CAAGTGAGAGAG right palm 2009     3  17 subject-1
#> 13    L3S360     CATCGTATCAAC right palm 2009     4  14 subject-1
#> 14    L5S104     CAGTGTCAGGAC     tongue 2008    10  28 subject-1
#> 15    L5S155     ATCTTAGACTGC     tongue 2009     1  20 subject-1
#> 16    L5S174     CAGACATTGCGT     tongue 2009     2  17 subject-1
#> 17    L5S203     CGATGCACCAGA     tongue 2009     3  17 subject-1
#> 18    L5S222     CTAGAGACTCTT     tongue 2009     4  14 subject-1
#> 19    L1S140     ATGGCAGCTCTA        gut 2008    10  28 subject-2
#> 20    L1S208     CTGAGATACGCG        gut 2009     1  20 subject-2
#> 21    L1S257     CCGACTGAGATG        gut 2009     3  17 subject-2
#> 22    L1S281     CCTCTCGTGATC        gut 2009     4  14 subject-2
#> 23    L2S240     CATATCGCAGTT  left palm 2008    10  28 subject-2
#> 24    L2S309     CGTGCATTATCA  left palm 2009     1  20 subject-2
#> 25    L2S357     CTAACGCAGTCA  left palm 2009     3  17 subject-2
#> 26    L2S382     CTCAATGACTCA  left palm 2009     4  14 subject-2
#> 27    L3S378     ATCGATCTGTGG right palm 2008    10  28 subject-2
#> 28     L4S63     CTCGTGGAGTAG right palm 2009     1  20 subject-2
#> 29    L4S112     GCGTTACACACA right palm 2009     3  17 subject-2
#> 30    L4S137     GAACTGTATCTC right palm 2009     4  14 subject-2
#> 31    L5S240     CTGGACTCATAG     tongue 2008    10  28 subject-2
#> 32     L6S20     GAGGCTCATCAT     tongue 2009     1  20 subject-2
#> 33     L6S68     GATACGTCCTGA     tongue 2009     3  17 subject-2
#> 34     L6S93     GATTAGCACTCT     tongue 2009     4  14 subject-2
#>    reported-antibiotic-usage days-since-experiment-start
#> 1                        Yes                           0
#> 2                         No                          84
#> 3                         No                         112
#> 4                         No                         140
#> 5                         No                          84
#> 6                         No                         112
#> 7                         No                         140
#> 8                         No                         168
#> 9                        Yes                           0
#> 10                        No                          84
#> 11                        No                         112
#> 12                        No                         140
#> 13                        No                         168
#> 14                       Yes                           0
#> 15                        No                          84
#> 16                        No                         112
#> 17                        No                         140
#> 18                        No                         168
#> 19                       Yes                           0
#> 20                        No                          84
#> 21                        No                         140
#> 22                        No                         168
#> 23                       Yes                           0
#> 24                        No                          84
#> 25                        No                         140
#> 26                        No                         168
#> 27                       Yes                           0
#> 28                        No                          84
#> 29                        No                         140
#> 30                        No                         168
#> 31                       Yes                           0
#> 32                        No                          84
#> 33                        No                         140
#> 34                        No                         168


## ------------------------------------------------
## Method `dataset$get_sample_tree`
## ------------------------------------------------


 tree <- ape::read.tree(strollur_example(
  "final.opti_mcc.jclass.ave.tre"))

 df <- read_mothur_shared(strollur_example("final.opti_mcc.shared"))

 data <- dataset$new("my_dataset")

 # assign abundance 'otu' bins
 assign(data = data, table = df, type = "bins", bin_type = "otu")
#> ℹ Assigned 531 otu bins.
#> [1] 531

 data$add_sample_tree(tree)
 data$get_sample_tree()
#> 
#> Phylogenetic tree with 19 tips and 18 internal nodes.
#> 
#> Tip labels:
#>   F3D9, F3D8, F3D6, F3D5, F3D2, F3D1, ...
#> 
#> Rooted; includes branch length(s).


## ------------------------------------------------
## Method `dataset$get_sequence_tree`
## ------------------------------------------------


 data <- dataset$new("my_dataset")
 tree <- ape::read.tree(strollur_example("final.phylip.tre.gz"))
 data$add_sequence_tree(tree)
#> ℹ Added 2425 sequences.
 data$get_sequence_tree()
#> 
#> Phylogenetic tree with 2425 tips and 2424 internal nodes.
#> 
#> Tip labels:
#>   M00967_43_000000000-A3JHG_1_1114_15727_25995, M00967_43_000000000-A3JHG_1_2109_19976_22044, M00967_43_000000000-A3JHG_1_1102_9244_9305, M00967_43_000000000-A3JHG_1_2101_14159_9619, M00967_43_000000000-A3JHG_1_1111_12315_7486, M00967_43_000000000-A3JHG_1_1107_12586_26826, ...
#> 
#> Rooted; includes branch length(s).
```
