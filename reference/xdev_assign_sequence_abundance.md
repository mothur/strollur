# xdev_assign_sequence_abundance

Assign sequence abundance and optionally assign sample and treatment
data to a
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Usage

``` r
xdev_assign_sequence_abundance(
  data,
  table,
  sequence_name = "sequence_name",
  abundance = "abundance",
  sample = "sample",
  treatment = "treatment",
  verbose = TRUE
)
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object

- table, :

  a data.frame containing sequence abundance assignments

- sequence_name, :

  a string containing the name of the column in 'table' that contains
  the sequence names. Default column name is 'sequence_name'.

- abundance, :

  a string containing the name of the column in 'table' that contains
  the sequence abundances. Default column name is 'abundance'.

- sample, :

  a string containing the name of the column in 'table' that contains
  the sequence samples. Default column name is 'sample'. (Optional)

- treatment, :

  a string containing the name of the column in 'table' that contains
  the sequence treatments. Default column name is 'treatment'.

- verbose, :

  a boolean whether or not you want progress messages. Default = TRUE.

## Value

an updated
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Examples

``` r

data <- new_dataset("my_dataset")
sequence_abundance <- readRDS(strollur_example("miseq_abundance_by_sample.rds"))

xdev_assign_sequence_abundance(data = data, table = sequence_abundance)
#> Assigned 2425 sequence abundances.
#> my_dataset:
#> 
#> 
#> Number of unique seqs: 2425 
#> Total number of seqs: 113963 
#> 
#> Total number of samples: 19 
#> Total number of treatments: 2 
#> 
```
