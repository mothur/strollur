# Merge bins in your [strollur object](https://mothur.org/strollur/reference/strollur.html)

Designed with package integration in mind, the merge bins function
allows you to merge bins in a [strollur
object](https://mothur.org/strollur/reference/strollur.html)

## Usage

``` r
xdev_merge_bins(data, bin_names, reason = "merged", bin_type = "otu")
```

## Arguments

- data:

  a [strollur
  object](https://mothur.org/strollur/reference/strollur.html).

- bin_names:

  a vector of strings containing the names of the bins you would like
  merge. The resulting merged bin will be stored in the first bin_id in
  the vector.

- reason:

  a string indicating why you are merging bins. Default = "merged".

- bin_type:

  a string indicating the type of bin clusters. Default = "otu"

## Value

an updated [strollur
object](https://mothur.org/strollur/reference/strollur.html)

## Examples

``` r

 data <- strollur::miseq_sop_example()
#> Added 2425 sequences.
#> Assigned 2425 sequence abundances.
#> Assigned 2425 sequence taxonomies.
#> Assigned 531 otu bins.
#> Assigned 2425 asv bins.
#> Assigned 63 phylotype bins.
#> Assigned 19 samples to treatments.
#> Assigned 171 samples distances.
#> Assigned 531 otu bin taxonomies.
#> Assigned 531 otu bin representative sequences.
#> Added a metadata report.
#> Added 2 resource references.
#> Added a contigs_report report.

 # to merge otu5 and otu6

 bins_to_merge <- c("Otu005", "Otu006")

 strollur::xdev_merge_bins(data,
                             bin_names = bins_to_merge)
#> miseq_sop:
#> 
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  375    249      0        3     0      1.00
#> 2.5%-tile:       1  375    252      0        4     0   2850.05
#> 25%-tile:        1  375    252      0        4     0  28491.50
#> Median:          1  375    253      0        4     0  56982.00
#> 75%-tile:        1  375    253      0        5     0  85472.50
#> 97.5%-tile:      1  375    254      0        6     0 111113.95
#> Maximum:         1  375    256      0        6     0 113963.00
#> Mean:            1  375    252      0        4     0  56982.00
#> 
#> scrap_summary:
#>   type trash_code unique total
#> 1  otu     merged      1  6621
#> 
#> Number of unique seqs: 2425 
#> Total number of seqs: 113963 
#> 
#> Total number of samples: 19 
#> Total number of treatments: 2 
#> Total number of otus: 530 
#> Total number of otu bin classifications: 530 
#> Total number of asvs: 2425 
#> Total number of asv bin classifications: 2425 
#> Total number of phylotypes: 63 
#> Total number of phylotype bin classifications: 63 
#> Total number of sequence classifications: 2425 
#> Total number of resource references: 2 
#> Total number of custom reports: 2 
#> 

 # If you look at the scrap report, you will see Otu006 with the trash code
 # set to "merged".

 strollur::report(data, type = "bin_scrap")
#>   bin_name trash_code
#> 1   Otu006     merged
```
