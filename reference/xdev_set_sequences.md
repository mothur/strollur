# xdev_set_sequences

Designed with package integration in mind, the set sequences function
allows you to change the nucleotide strings of sequences in a
[strollur](https://mothur.org/strollur/reference/strollur.html) object.
For example, set_sequences may be used after alignment to overwrite the
unaligned sequences with aligned sequences.

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

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object

- sequence_names, :

  a vector of strings containing sequence names

- sequences, :

  a vector of strings containing sequence nucleotide strings

- comments, :

  a vector of strings containing sequence comments. (Optional)

## Value

an updated
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Examples

``` r

data <- new_dataset(dataset_name = "my_dataset")

xdev_add_sequences(data = data,
              table = data.frame(sequence_name = c("seq1", "seq2",
                                                  "seq3", "seq4")))
#> Added 4 sequences.
#> my_dataset:
#> 
#> 
#> Number of unique seqs: 4 
#> Total number of seqs: 4 
#> 
#> 

xdev_set_sequences(data = data,
                   sequence_names = c("seq1", "seq2","seq3", "seq4"),
                   sequences = c("ATTGC", "ACTGC", "AGTGC", "TTTGC"))
#> my_dataset:
#> 
#>             starts ends nbases ambigs polymers numns numseqs
#> Minimum:         1    5      5      0        1     0    1.00
#> 2.5%-tile:       1    5      5      0        1     0    0.10
#> 25%-tile:        1    5      5      0        1     0    1.00
#> Median:          1    5      5      0        1     0    2.00
#> 75%-tile:        1    5      5      0        2     0    3.00
#> 97.5%-tile:      1    5      5      0        2     0    3.90
#> Maximum:         1    5      5      0        3     0    4.00
#> Mean:            1    5      5      0        1     0    2.14
#> 
#> Number of unique seqs: 4 
#> Total number of seqs: 4 
#> 
#> 
```
