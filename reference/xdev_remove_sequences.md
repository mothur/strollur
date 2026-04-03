# xdev_remove_sequences

Designed with package integration in mind, the remove sequences function
allows you to remove sequences from a
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Usage

``` r
xdev_remove_sequences(data, sequence_names, trash_tags)
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object.

- sequence_names, :

  vector of strings containing the names of the sequences to remove

- trash_tags:

  vector of strings containing the reasons for the sequences removals

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

# For the sake of example let's remove the first 3 sequences from
# miseq_sop_example:

seqs_to_remove <- c("M00967_43_000000000-A3JHG_1_2101_16474_12783",
                   "M00967_43_000000000-A3JHG_1_1113_12711_3318",
                   "M00967_43_000000000-A3JHG_1_2108_14707_9807")
trash_codes <- c("example", "removing", "sequences")

xdev_remove_sequences(data = data, sequence_names = seqs_to_remove,
                      trash_tags = trash_codes)

# If you look at the scrap report, you the sequences names, listed with the
# trash codes set to "example", "removing", "sequences".

report(data = data, type = "sequence_scrap")
#>                                             id trash_code
#> 1 M00967_43_000000000-A3JHG_1_2101_16474_12783    example
#> 2  M00967_43_000000000-A3JHG_1_1113_12711_3318   removing
#> 3  M00967_43_000000000-A3JHG_1_2108_14707_9807  sequences

# You can see from the get_num_sequences function that the removed
# sequence's abundances are removed from the dataset.

count(data = data, type = "sequences")
#> [1] 113960
```
