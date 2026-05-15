# xdev_assign_treatments

Assign samples to treatments in a
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Usage

``` r
xdev_assign_treatments(
  data,
  table,
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

  a data.frame containing sample treatment assignments

- sample, :

  a string containing the name of the column in 'table' that contains
  the samples. Default column name is 'sample'.

- treatment, :

  a string containing the name of the column in 'table' that contains
  the treatments. Default column name is 'treatment'.

- verbose, :

  a boolean indicating whether or not you want progress messages.
  Default = TRUE.

## Value

an updated
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Examples

``` r

data <- new_dataset("my_dataset")
sequence_abundance <- readRDS(strollur_example("miseq_abundance_by_sample.rds"))

xdev_assign_sequence_abundance(data, sequence_abundance)
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

sample_assignments <- readRDS(strollur_example("miseq_sample_design.rds"))

xdev_assign_treatments(data, sample_assignments)
#> Assigned 19 samples to treatments.
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
