# xdev_set_sequences

Designed with package integration in mind, the set sequences function
allows you to change the nucleotide strings of sequences in a
[dataset](dataset.md) object. For example, set_sequences may be used
after alignment to overwrite the unaligned sequences with aligned
sequences.

## Usage

``` r
xdev_set_sequences(
  data,
  sequence_names,
  sequences,
  comments = as.character(c())
)
```

## Arguments

- data, :

  a [dataset](dataset.md) object

- sequence_names, :

  a vector of strings containing sequence names

- sequences, :

  a vector of strings containing sequence nucleotide strings

- comments, :

  a vector of strings containing sequence comments. (Optional)

## Examples

``` r
data <- new_dataset(dataset_name = "my_dataset")

xdev_add_sequences(data = data,
              table = data.frame(sequence_names = c("seq1", "seq2",
                                                  "seq3", "seq4")))
#> ℹ Added 4 sequences.
#> [1] 4

xdev_set_sequences(data = data,
                   sequence_names = c("seq1", "seq2","seq3", "seq4"),
                   sequences = c("ATTGC", "ACTGC", "AGTGC", "TTTGC"))
```
