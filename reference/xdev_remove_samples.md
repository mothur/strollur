# xdev_remove_samples

Designed with package integration in mind, the remove samples function
allows you to remove samples from a
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Usage

``` r
xdev_remove_samples(data, samples, reason = "remove_samples")
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object.

- samples, :

  vector of strings containing the names of the samples to remove.

- reason, :

  string containing the reason for removal. Default = "remove_samples".

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

count(data = data, type = "sample")
#> [1] 19

# To remove samples 'F3D0' and 'F3D1'

xdev_remove_samples(data, c("F3D0", "F3D1"))
#> miseq_sop:
#> 
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  375    249      0        3     0      1.00
#> 2.5%-tile:       1  375    252      0        4     0   2579.00
#> 25%-tile:        1  375    252      0        4     0  25781.00
#> Median:          1  375    253      0        4     0  51561.00
#> 75%-tile:        1  375    253      0        5     0  77341.00
#> 97.5%-tile:      1  375    254      0        6     0 100543.00
#> Maximum:         1  375    256      0        6     0 103121.00
#> Mean:            1  375    252      0        4     0  51561.00
#> 
#> scrap_summary:
#>        type     trash_code unique total
#> 1  sequence remove_samples    204   215
#> 2       otu remove_samples     30    31
#> 3       asv remove_samples    204   215
#> 4 phylotype remove_samples      4     4
#> 
#> Number of unique seqs: 2221 
#> Total number of seqs: 103120 
#> 
#> Total number of samples: 17 
#> Total number of treatments: 2 
#> Total number of otus: 501 
#> Total number of otu bin classifications: 501 
#> Total number of asvs: 2221 
#> Total number of asv bin classifications: 2221 
#> Total number of phylotypes: 59 
#> Total number of phylotype bin classifications: 59 
#> Total number of sequence classifications: 2221 
#> Total number of resource references: 2 
#> Total number of custom reports: 1 
#> Your dataset includes metadata 
#> 

count(data = data, type = "sample")
#> [1] 17
```
