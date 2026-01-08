# xdev_add_sequences

Add sequence data to a [dataset](dataset.md) object

## Usage

``` r
xdev_add_sequences(
  data,
  table,
  reference = NULL,
  sequence_name = "sequence_names",
  sequence = "sequences",
  comment = "comments",
  verbose = TRUE
)
```

## Arguments

- data, :

  a [dataset](dataset.md) object

- table, :

  a data.frame containing names, sequences(optional) and
  comments(optional).

- reference, :

  a list created by the function \[new_reference\]. Optional.

- sequence_name, :

  a string containing the name of the column in 'table' that contains
  the sequence names. Default column name is 'sequence_names'.

- sequence, :

  a string containing the name of the column in 'table' that contains
  the sequence nucleotide strings. Default column name is 'sequences'.

- comment, :

  a string containing the name of the column in 'table' that contains
  the sequence comments. Default column name is 'comments'.

- verbose, :

  a boolean whether or not you want progress messages. Default = TRUE.

## Value

double containing the number of sequences added

## Examples

``` r
 data <- new_dataset("miseq_sop", 2)
 fasta_data <- read_fasta(rdataset_example("final.fasta"))
 xdev_add_sequences(data, fasta_data)
#> ℹ Added 2425 sequences.
#> [1] 2425

# With the additional parameters to add information about the reference

 data <- new_dataset("miseq_sop", 2)
 fasta_data <- read_fasta(rdataset_example("final.fasta"))

 xdev_add_sequences(data, fasta_data,
               new_reference("silva.bacteria.fasta",
               "1.38.1",
               "alignment by mothur2 v1.0 using default options",
               "https://mothur.org/wiki/silva_reference_files/"))
#> ℹ Added 2425 sequences.
#> [1] 2425

# You can also add references using the 'add_references' function.
```
