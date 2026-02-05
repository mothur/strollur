# xdev_set_abundance

Designed with package integration in mind, the set abundance function
allows you to change the abundances of sequences in a
[dataset](https://mothur.org/strollur/reference/dataset.md) object
without samples.

## Usage

``` r
xdev_set_abundance(
  data,
  sequence_names,
  sequence_abundances,
  reason = "update"
)
```

## Arguments

- data, :

  a [dataset](https://mothur.org/strollur/reference/dataset.md) object

- sequence_names, :

  a vector of strings containing sequence names

- sequence_abundances, :

  vector containing the abundances of each sequence.

- reason, :

  a string containing the trash tag to be applied to any sequences set
  to 0 abundance. Default = "update".

## Examples

``` r
names <- c("seq1", "seq2", "seq3",  "seq4")
abunds <- c(1250, 65, 50, 4)

data <- new_dataset(dataset_name = "my_dataset")

xdev_assign_sequence_abundance(data = data, table = data.frame(sequence_names = names,
                                           abundances = abunds))
#> ℹ Assigned 4 sequence abundances.
#> [1] 4
abundance(data = data, type = "sequences")
#>   sequence_names abundances
#> 1           seq1       1250
#> 2           seq2         65
#> 3           seq3         50
#> 4           seq4          4

seqs_to_update <- c("seq1", "seq3")
new_abunds <- c(1000, 100)

xdev_set_abundance(data = data,
                   sequence_names = seqs_to_update,
                   sequence_abundances = new_abunds)

abundance(data = data, type = "sequences")
#>   sequence_names abundances
#> 1           seq1       1000
#> 2           seq2         65
#> 3           seq3        100
#> 4           seq4          4
```
