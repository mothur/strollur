# xdev_merge_bins

Designed with package integration in mind, the merge bins function
allows you to merge bins in a
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Usage

``` r
xdev_merge_bins(data, bin_names, reason = "merged", bin_type = "otu")
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object.

- bin_names, :

  a vector of strings containing the names of the bins you would like
  merge. The resulting merged bin will be stored in the first bin_id in
  the vector.

- reason, :

  a string indicating why you are merging bins. Default = "merged".

- bin_type, :

  a string indicating the type of bin clusters. Default = "otu"

## Value

an updated
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Examples

``` r

 data <- miseq_sop_example()
#> Added 2425 sequences.
#> Assigned 2425 sequence abundances.
#> Assigned 2425 sequence taxonomies.
#> Assigned 531 otu bins.
#> Assigned 2425 asv bins.
#> Assigned 63 phylotype bins.
#> Assigned 19 samples to treatments.
#> Assigned 531 otu bin taxonomies.
#> Assigned 531 otu bin representative sequences.
#> Added metadata.
#> Added 2 resource references.
#> Added a contigs_report.

 # to merge otu5 and otu6

 bins_to_merge <- c("Otu005", "Otu006")

 xdev_merge_bins(data = data, bin_names = bins_to_merge)
#> miseq_sop:
#> 
#>             starts ends   nbases ambigs polymers numns    numseqs
#> Minimum:         1  375 249.0000      0 3.000000     0      1.000
#> 2.5%-tile:       1  375 252.0000      0 4.000000     0   2850.075
#> 25%-tile:        1  375 252.0000      0 4.000000     0  28491.750
#> Median:          1  375 253.0000      0 4.000000     0  56982.500
#> 75%-tile:        1  375 253.0000      0 5.000000     0  85473.250
#> 97.5%-tile:      1  375 254.0000      0 6.000000     0 111114.925
#> Maximum:         1  375 256.0000      0 6.000000     0 113964.000
#> Mean:            1  375 252.7406      0 4.496082     0  56982.500
#>     Length Overlap_Length Overlap_Start Overlap_End MisMatches Num_Ns
#> 1 250.0000       232.0000      0.000000    248.0000   0.000000      0
#> 2 252.0000       246.0000      1.000000    250.0000   0.000000      0
#> 3 252.0000       249.0000      2.000000    251.0000   0.000000      0
#> 4 253.0000       249.0000      2.000000    251.0000   1.000000      0
#> 5 253.0000       250.0000      2.000000    251.0000   5.000000      0
#> 6 254.0000       251.0000      4.000000    253.0000  26.000000      0
#> 7 270.0000       255.0000     22.000000    256.0000 120.000000      0
#> 8 252.7575       249.1501      2.005361    251.1555   5.162474      0
#>   Expected_Errors
#> 1      1.00000000
#> 2      1.00000000
#> 3      1.00000000
#> 4      1.00000000
#> 5      1.00000000
#> 6      1.00000000
#> 7      4.00000000
#> 8      0.07385095
#>   type trash_code unique total
#> 1  otu     merged      1  6621
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  375    249      0        3     0      1.00
#> 2.5%-tile:       1  375    252      0        4     0   2850.08
#> 25%-tile:        1  375    252      0        4     0  28491.75
#> Median:          1  375    253      0        4     0  56982.50
#> 75%-tile:        1  375    253      0        5     0  85473.25
#> 97.5%-tile:      1  375    254      0        6     0 111114.93
#> Maximum:         1  375    256      0        6     0 113964.00
#> Mean:            1  375    252      0        4     0  56982.50
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
#> Total number of custom reports: 1 
#> Your dataset includes metadata 
#> 

 # If you look at the scrap report, you will see Otu006 with the trash code
 # set to "merged".

 report(data = data, type = "bin_scrap")
#>       id trash_code
#> 1 Otu006     merged
```
