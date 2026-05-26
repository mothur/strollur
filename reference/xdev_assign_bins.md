# xdev_assign_bins

Add bin assignments to a
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Usage

``` r
xdev_assign_bins(
  data,
  table,
  bin_type = "otu",
  reference = NULL,
  bin_name = "bin_name",
  abundance = "abundance",
  sample = "sample",
  sequence_name = "sequence_name",
  verbose = TRUE
)
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object

- table, :

  a data.frame containing bin_data assignments

- bin_type:

  a string indicating the type of bin assignments. Default "otu".

- reference, :

  a list created by the function \[new_reference\]. Optional.

- bin_name, :

  a string containing the name of the column in 'table' that contains
  the bin names. Default column name is 'bin_name'.

- abundance, :

  a string containing the name of the column in 'table' that contains
  the bin abundances. Default column name is 'abundance'. Note: You must
  provide either abundance or sequence_name in the table.

- sample, :

  a string containing the name of the column in 'table' that contains
  the sample names for datasets where the abundances are broken down by
  sample. Default column name is 'sample'.

- sequence_name, :

  a string containing the name of the column in 'table' that contains
  the sequence names. Default column name is 'sequence_name'. Note: You
  must provide either abundance or sequence_name in the table.

- verbose, :

  a boolean whether or not you want progress messages. Default = TRUE.

## Value

an updated
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Examples

``` r

  # To assign sequences to bins:

  data <- new_dataset(dataset_name = "miseq_sop")
  otu_data <- read_mothur_list(list = strollur_example("final.opti_mcc.list.gz"))

  xdev_assign_bins(data = data, table = otu_data, bin_type = "otu")
#> Assigned 531 otu bins.
#> miseq_sop:
#> 
#> data frame with 0 columns and 0 rows
#> 
#> Number of unique seqs: 2425 
#> Total number of seqs: 2425 
#> 
#> Total number of otus: 531 
#> 

  # To add abundance only bin assignments:

  data <- new_dataset(dataset_name = "miseq_sop")
  otu_data <- read_mothur_rabund(rabund = strollur_example("final.opti_mcc.rabund"))

  xdev_assign_bins(data = data, table = otu_data, bin_type = "otu")
#> Assigned 531 otu bins.
#> miseq_sop:
#> 
#> data frame with 0 columns and 0 rows
#> 
#> Number of unique seqs: 531 
#> Total number of seqs: 2425 
#> 
#> Total number of otus: 531 
#> 

  # To add abundance bin assignments parsed by sample:

  data <- new_dataset(dataset_name = "miseq_sop")
  otu_data <- readRDS(strollur_example("miseq_shared_otu.rds"))

  xdev_assign_bins(data = data, table = otu_data, bin_type = "otu")
#> Assigned 531 otu bins.
#> miseq_sop:
#> 
#> data frame with 0 columns and 0 rows
#> 
#> Number of unique seqs: 531 
#> Total number of seqs: 113963 
#> 
#> Total number of samples: 19 
#> Total number of otus: 531 
#> 
```
