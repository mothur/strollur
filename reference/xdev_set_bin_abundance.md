# xdev_set_bin_abundance

Designed with package integration in mind, the set bin abundance
function allows you to change the abundances of bins in a
[dataset](dataset.md) object without sample data.

## Usage

``` r
xdev_set_bin_abundance(
  data,
  bin_names,
  abundances,
  type = "otu",
  reason = "update"
)
```

## Arguments

- data, :

  a [dataset](dataset.md) object

- bin_names, :

  a vector strings containing of bin names to set the abundances for.

- abundances, :

  vector containing the abundances of each bin.

- type, :

  a string indicating the type of clusters. Default = "otu".

- reason, :

  a string containing the trash tag to be applied to any bins set to 0
  abundance. Default = "update".

## Examples

``` r
  # For example sake, let's create a dataset with 3 bins:

  data <- new_dataset(dataset_name = "my_dataset")

  bin_ids <- c("bin1", "bin2", "bin3")
  abundances <- c(110, 525, 80)

  xdev_assign_bins(data = data, table = data.frame(bin_names = bin_ids,
                                              abundances = abundances))
#> ℹ Assigned 3 otu bins.
#> [1] 3

  abundance(data, type = "bins")
#>   otu_id abundance
#> 1   bin1       110
#> 2   bin2       525
#> 3   bin3        80

  # Now we can use set_bin_abundance to change the abundances of bin1 and
  # bin2

  bins <- c("bin1", "bin2")
  new_abunds <- c(300, 250)

  xdev_set_bin_abundance(data = data,
                         bin_names = bins,
                         abundances = new_abunds)

  abundance(data, type = "bins")
#>   otu_id abundance
#> 1   bin1       300
#> 2   bin2       250
#> 3   bin3        80
```
