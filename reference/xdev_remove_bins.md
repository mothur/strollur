# xdev_remove_bins

Designed with package integration in mind, the remove bins function
allows you to remove bins from a
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Usage

``` r
xdev_remove_bins(data, bin_names, trash_tags, bin_type = "otu")
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object.

- bin_names, :

  a vector of strings containing the names of the bins you would like
  removed.

- trash_tags, :

  a vector of strings containing the reasons you are removing each bin

- bin_type:

  a string indicating the type of clusters.

## Value

an updated
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Examples

``` r

  data <- new_dataset(dataset_name = "my_dataset")

  bin_names <- c("bin1", "bin2", "bin3")
  abundances <- c(110, 525, 80)

  xdev_assign_bins(data = data,
                   table = data.frame(bin_name = bin_names,
                                      abundance = abundances),
                   bin_type = "otu")
#> Assigned 3 otu bins.
#> my_dataset:
#> 
#> 
#> Number of unique seqs: 3 
#> Total number of seqs: 715 
#> 
#> Total number of otus: 3 
#> 

  count(data = data, type = "bin", bin_type = "otu")
#> [1] 3

  bins_to_remove <- c("bin1")
  trash_tag <- c("bad_bin")

  xdev_remove_bins(data = data,
                   bin_names = bins_to_remove,
                   trash_tags = trash_tag)
#> my_dataset:
#> 
#> scrap_summary:
#>       type trash_code unique total
#> 1 sequence    bad_bin      1   110
#> 2      otu    bad_bin      1   110
#> 
#> Number of unique seqs: 2 
#> Total number of seqs: 605 
#> 
#> Total number of otus: 2 
#> 

  count(data = data, type = "bin", bin_type = "otu")
#> [1] 2
```
