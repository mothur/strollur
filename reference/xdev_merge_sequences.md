# xdev_merge_sequences

Designed with package integration in mind, the merge sequences function
allows you to merge sequences in a
[strollur](https://mothur.org/strollur/reference/strollur.md) object.

## Usage

``` r
xdev_merge_sequences(data, sequence_names, reason = "merged")
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.md)
  object.

- sequence_names, :

  a vector of strings containing the names of the sequences you would
  like merge. The resulting merged sequence will be stored in the first
  sequence name in the vector.

- reason:

  a string indicating why you are merging sequences. Default = "merged"

## Examples

``` r
sequence_names <- c("seq1", "seq2", "seq3", "seq3",
               "seq4", "seq4", "seq5", "seq6",
               "seq7", "seq8", "seq9", "seq9",
               "seq10", "seq10", "seq10", "seq10")

samples <- c("sample1", "sample2", "sample4", "sample5",
             "sample1", "sample2", "sample1", "sample1",
             "sample2", "sample4", "sample4", "sample5",
             "sample1", "sample3", "sample5", "sample6")

abundances <- c(10, 10, 5, 5, 5, 5,
                10, 10, 10, 10, 5, 5,
                1, 2, 3, 4)

data <- new_dataset("my_data")


assign(data = data,
       table = data.frame(sequence_names = sequence_names,
                          abundances = abundances,
                          samples = samples),
       type = "sequence_abundance")
#> ℹ Assigned 10 sequence abundances.
#> [1] 10

# For the sake of example let's merge the first 3 sequences.

seqs_to_merge <- c("seq1", "seq2", "seq3")

xdev_merge_sequences(data = data, sequence_names = seqs_to_merge)

# If you look at the scrap report, you will see the second two sequence
# names, listed with the trash code set to "merged".

report(data = data, type = "sequence_scrap")
#>     id trash_code
#> 1 seq2     merged
#> 2 seq3     merged

# You can see from the get_num_sequences function that the merged sequence's
# abundances are added to the first sequence.

count(data = data, type = "sequences")
#> [1] 100
```
