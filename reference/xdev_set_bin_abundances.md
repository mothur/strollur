# xdev_set_bin_abundances

Designed with package integration in mind, the set bin abundances
function allows you to change the abundances of bins in a
[dataset](dataset.md) object with sample data.

## Usage

``` r
xdev_set_bin_abundances(
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

  2D vector (\[num_seqs\]\[num_samples\]) containing the abundances of
  each bin parsed by sample.

- type:

  a string indicating the type of clusters. Default = "otu".

- reason, :

  a string containing the trash tag to be applied to any bins set to 0
  abundance. Default = "update".

## Examples

``` r
  # For example sake, let's create a dataset with 3 bins:

  data <- new_dataset(dataset_name = "my_dataset")

  bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
  samples <- c("sample1", "sample2", "sample5", "sample1", "sample3",
               "sample1")
  sample_abundances <- c(10, 100, 1, 500, 25, 80)

  xdev_assign_bins(data = data, table = data.frame(bin_names = bin_ids,
                                              abundances = sample_abundances,
                                              samples = samples))
#> ℹ Assigned 3 otu bins.
#> [1] 3

  # You can see bin1's abundances parsed by sample using abundance:
  abundance(data, type = "bins", by_sample = TRUE)
#>   bin_names abundances samples
#> 1      bin1         10 sample1
#> 2      bin1        100 sample2
#> 3      bin1          1 sample5
#> 4      bin2        500 sample1
#> 5      bin2         25 sample3
#> 6      bin3         80 sample1

  # You can change bin1's abundances as follows:

  new_bin1_abunds <- list(c(10,50,0,0))
  bins <- c("bin1")

  xdev_set_bin_abundances(data = data,
                          bin_names = bins,
                          abundances = new_bin1_abunds)

  abundance(data, type = "bins", by_sample = TRUE)
#>   bin_names abundances samples
#> 1      bin1         10 sample1
#> 2      bin1         50 sample2
#> 3      bin2        500 sample1
#> 4      bin2         25 sample3
#> 5      bin3         80 sample1
```
