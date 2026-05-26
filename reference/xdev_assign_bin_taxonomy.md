# xdev_assign_bin_taxonomy

Assign bin classifications to a
[strollur](https://mothur.org/strollur/reference/strollur.html) object

Note, if you assign sequence taxonomies and assign bins, 'Dataset' will
find the concensus taxonomy for each bin for you.

## Usage

``` r
xdev_assign_bin_taxonomy(
  data,
  table,
  bin_type = "otu",
  reference = NULL,
  bin_name = "bin_name",
  taxonomy = "taxonomy",
  verbose = TRUE
)
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object

- table, :

  a data.frame containing bin taxonomy assignments

- bin_type:

  a string indicating the type of bin assignments. Default "otu".

- reference, :

  a list created by the function \[new_reference\]. Optional.

- bin_name, :

  a string containing the name of the column in 'table' that contains
  the bin names. Default column name is 'bin_name'.

- taxonomy, :

  a string containing the name of the column in 'table' that contains
  the bin taxonomies. Default column name is 'taxonomy'.

- verbose, :

  a boolean whether or not you want progress messages. Default = TRUE.

## Value

an updated
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Examples

``` r

otu_data <- read_mothur_cons_taxonomy(strollur_example(
                        "final.cons.taxonomy"))

data <- new_dataset(dataset_name = "my_dataset")

# assign otu abundances
xdev_assign_bins(data = data, table = otu_data, bin_type = "otu")
#> Assigned 531 otu bins.
#> my_dataset:
#> 
#> data frame with 0 columns and 0 rows
#> 
#> Number of unique seqs: 531 
#> Total number of seqs: 113963 
#> 
#> Total number of otus: 531 
#> 

# assign otu classifications
xdev_assign_bin_taxonomy(data = data, table = otu_data,
                         bin_type = "otu")
#> Assigned 531 otu bin taxonomies.
#> my_dataset:
#> 
#> data frame with 0 columns and 0 rows
#> 
#> Number of unique seqs: 531 
#> Total number of seqs: 113963 
#> 
#> Total number of otus: 531 
#> Total number of otu bin classifications: 531 
#> 
```
