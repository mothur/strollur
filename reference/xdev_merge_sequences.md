# xdev_merge_sequences

Designed with package integration in mind, the merge sequences function
allows you to merge sequences in a [dataset](dataset.md) object.

## Usage

``` r
xdev_merge_sequences(data, sequence_names, reason = "merged")
```

## Arguments

- data, :

  a [dataset](dataset.md) object.

- sequence_names, :

  a vector of strings containing the names of the sequences you would
  like merge. The resulting merged sequence will be stored in the first
  sequence name in the vector.

- reason:

  a string indicating why you are merging sequences. Default = "merged"

## Examples

``` r
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

count(data = data, type = "sequences")
#> [1] 113963

# For the sake of example let's merge the first 3 sequences from
# miseq_sop_example:

seqs_to_merge <- c("M00967_43_000000000-A3JHG_1_2101_16474_12783",
                   "M00967_43_000000000-A3JHG_1_1113_12711_3318",
                   "M00967_43_000000000-A3JHG_1_2108_14707_9807")

xdev_merge_sequences(data = data, sequence_names = seqs_to_merge)

# If you look at the scrap report, you will see the second two sequence
# names, listed with the trash code set to "merged".

report(data = data, type = "sequence_scrap")
#>                                            id trash_code
#> 1 M00967_43_000000000-A3JHG_1_1113_12711_3318     merged
#> 2 M00967_43_000000000-A3JHG_1_2108_14707_9807     merged

# You can see from the get_num_sequences function that the merged sequence's
# abundances are added to the first sequence.

count(data = data, type = "sequences")
#> [1] 113963
```
