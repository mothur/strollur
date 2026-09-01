# Set abundances of sequences in a [strollur object](https://mothur.org/strollur/reference/strollur.html) with sample data

Designed with package integration in mind, the set abundances function
allows you to change the abundances of sequences in a [strollur
object](https://mothur.org/strollur/reference/strollur.html) with
samples.

## Usage

``` r
xdev_set_abundances(data, sequence_names, abundances, reason = "update")
```

## Arguments

- data:

  a [strollur
  object](https://mothur.org/strollur/reference/strollur.html)

- sequence_names:

  a vector of strings containing sequence names

- abundances:

  2D vector num_seqs x num_samples containing the abundances of each
  sequence parsed by sample.

- reason:

  a string containing the trash tag to be applied to any sequences set
  to 0 abundance. Default = "update".

## Value

an updated [strollur
object](https://mothur.org/strollur/reference/strollur.html)

## Examples

``` r

data <- strollur::new_dataset(dataset_name = "my_dataset")

sequence_names <- c("seq1", "seq1", "seq1", "seq2", "seq2", "seq2", "seq3",
                    "seq3", "seq4")
samples <- c("sample2", "sample3", "sample4", "sample2", "sample3",
             "sample4", "sample2", "sample3", "sample4")
abundances <- c(250, 400, 500, 25, 40, 50, 25, 25, 4)
table <- data.frame(sequence_name = sequence_names,
                    abundance = abundances,
                    sample = samples)

strollur::xdev_assign_sequence_abundance(data,
                                           table = table)
#> Assigned 4 sequence abundances.
#> my_dataset:
#> 
#> 
#> Number of unique seqs: 4 
#> Total number of seqs: 1319 
#> 
#> Total number of samples: 3 
#> 

seqs_to_update <- c("seq4")
new_abunds <- list(c(20, 10, 4))

strollur::xdev_set_abundances(data,
                    sequence_names = seqs_to_update,
                    abundances = new_abunds)
#> my_dataset:
#> 
#> 
#> Number of unique seqs: 4 
#> Total number of seqs: 1349 
#> 
#> Total number of samples: 3 
#> 
```
