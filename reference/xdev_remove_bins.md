# xdev_remove_bins

Designed with package integration in mind, the remove bins function
allows you to remove bins from a
[dataset](https://mothur.org/strollur/reference/dataset.md) object

## Usage

``` r
xdev_remove_bins(data, bin_names, trash_tags, bin_type = "otu")
```

## Arguments

- data, :

  a [dataset](https://mothur.org/strollur/reference/dataset.md) object.

- bin_names, :

  a vector of strings containing the names of the bins you would like
  removed.

- trash_tags, :

  a vector of strings containing the reasons you are removing each bin

- bin_type:

  a string indicating the type of clusters.

## Examples

``` r
  data <- new_dataset(dataset_name = "my_dataset")

  bin_names <- c("bin1", "bin2", "bin3")
  abundances <- c(110, 525, 80)

  xdev_assign_bins(data = data, table = data.frame(bin_names = bin_names,
                               abundances = abundances), bin_type = "otu")
#> ℹ Assigned 3 otu bins.
#> [1] 3

  count(data = data, type = "bins", bin_type = "otu")
#> [1] 3

  bins_to_remove <- c("bin1")
  trash_tag <- c("bad_bin")

  xdev_remove_bins(data = data,
                   bin_names = bins_to_remove,
                   trash_tags = trash_tag)

  count(data = data, type = "bins", bin_type = "otu")
#> [1] 2
```
