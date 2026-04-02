# xdev_merge_bins

Designed with package integration in mind, the merge bins function
allows you to merge bins in a
[strollur](https://mothur.org/strollur/reference/strollur.md) object

## Usage

``` r
xdev_merge_bins(data, bin_names, reason = "merged", bin_type = "otu")
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.md)
  object.

- bin_names, :

  a vector of strings containing the names of the bins you would like
  merge. The resulting merged bin will be stored in the first bin_id in
  the vector.

- reason, :

  a string indicating why you are merging bins. Default = "merged".

- bin_type, :

  a string indicating the type of bin clusters. Default = "otu"

## Examples

``` r
 data <- miseq_sop_example()
#> ℹ Added 2425 sequences.
#> ℹ Assigned 2425 sequence abundances.
#> ℹ Assigned 2425 sequence taxonomies.
#> ℹ Assigned 531 otu bins.
#> ℹ Assigned 2425 asv bins.
#> ℹ Assigned 63 phylotype bins.
#> ℹ Assigned 19 samples to treatments.
#> ℹ Assigned 531 otu bin taxonomies.
#> ℹ Assigned 531 otu bin representative sequences.
#> ℹ Added metadata.
#> ℹ Added 2 resource references.
#> ℹ Added a contigs_report.

 # to merge otu5 and otu6

 bins_to_merge <- c("Otu005", "Otu006")

 xdev_merge_bins(data = data, bin_names = bins_to_merge)

 # If you look at the scrap report, you will see Otu006 with the trash code
 # set to "merged".

 report(data = data, type = "bin_scrap")
#>       id trash_code
#> 1 Otu006     merged
```
