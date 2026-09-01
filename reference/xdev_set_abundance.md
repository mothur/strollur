# Set abundances of sequences in a [strollur object](https://mothur.org/strollur/reference/strollur.html) without sample data

Designed with package integration in mind, the set abundance function
allows you to change the abundances of sequences in a
[strollur](https://mothur.org/strollur/reference/strollur.html) object
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

- data:

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object

- sequence_names:

  a vector of strings containing sequence names

- sequence_abundances:

  vector containing the abundances of each sequence.

- reason:

  a string containing the trash tag to be applied to any sequences set
  to 0 abundance. Default = "update".

## Value

an updated
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Examples

``` r

names <- c("seq1", "seq2", "seq3",  "seq4")
abunds <- c(1250, 65, 50, 4)
table <- data.frame(sequence_name = names,
                   abundance = abunds)

data <- strollur::new_dataset(dataset_name = "my_dataset")

strollur::xdev_assign_sequence_abundance(data,
                                         table = table)
#> Assigned 4 sequence abundances.
#> my_dataset:
#> 
#> 
#> Number of unique seqs: 4 
#> Total number of seqs: 1369 
#> 
#> 
strollur::abundance(data, type = "sequence")
#>   sequence_name abundance
#> 1          seq1      1250
#> 2          seq2        65
#> 3          seq3        50
#> 4          seq4         4

seqs_to_update <- c("seq1", "seq3")
new_abunds <- c(1000, 100)

strollur::xdev_set_abundance(data,
                   sequence_names = seqs_to_update,
                   sequence_abundances = new_abunds)
#> my_dataset:
#> 
#> 
#> scrap_summary:
#>       type trash_code unique total
#> 1 sequence     update      0   200
#> 
#> Number of unique seqs: 4 
#> Total number of seqs: 1169 
#> 
#> 

strollur::abundance(data, type = "sequence")
#>   sequence_name abundance
#> 1          seq1      1000
#> 2          seq2        65
#> 3          seq3       100
#> 4          seq4         4
```
