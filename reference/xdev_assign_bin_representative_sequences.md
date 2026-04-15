# xdev_assign_bin_representative_sequences

Assign representative sequences to bins.

## Usage

``` r
xdev_assign_bin_representative_sequences(
  data,
  table,
  bin_type = "otu",
  reference = NULL,
  bin_name = "bin_names",
  sequence_name = "sequence_names",
  verbose = TRUE
)
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object

- table, :

  a data.frame containing bin representative assignments

- bin_type:

  a string indicating the type of bin assignments. Default "otu".

- reference, :

  a list created by the function \[new_reference\]. Optional.

- bin_name, :

  a string containing the name of the column in 'table' that contains
  the bin names. Default column name is 'bin_names'.

- sequence_name:

  a string containing the name of the column in 'table' that contains
  the bin names. Default column name is 'sequence_names'.

- verbose, :

  a boolean whether or not you want progress messages. Default = TRUE.

## Value

double containing the number of representative sequences assigned

## Examples

``` r
  miseq <- miseq_sop_example()
#> Added 2425 sequences.
#> Assigned 2425 sequence abundances.
#> Assigned 2425 sequence taxonomies.
#> Assigned 531 otu bins.
#> Assigned 2425 asv bins.
#> Assigned 63 phylotype bins.
#> Assigned 19 samples to treatments.
#> Assigned 531 otu bin taxonomies.
#> Assigned 531 otu bin representative sequences.
#> Added metadata.
#> Added 2 resource references.
#> Added a contigs_report.

  bin_reps <- readRDS(strollur_example(
         "miseq_representative_sequences.rds"))

  xdev_assign_bin_representative_sequences(data = miseq,
                                      table = bin_reps,
                                      bin_type = "otu")
#> Assigned 531 otu bin representative sequences.
#> [1] 531
```
