# xdev_remove_sequences

Designed with package integration in mind, the remove sequences function
allows you to remove sequences from a
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Usage

``` r
xdev_remove_sequences(data, sequence_names, trash_tags)
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object.

- sequence_names, :

  vector of strings containing the names of the sequences to remove

- trash_tags:

  vector of strings containing the reasons for the sequences removals

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

count(data = data, type = "sequence")
#> [1] 113963

# For the sake of example let's remove the first 3 sequences from
# miseq_sop_example:

seqs_to_remove <- c("M00967_43_000000000-A3JHG_1_2101_16474_12783",
                   "M00967_43_000000000-A3JHG_1_1113_12711_3318",
                   "M00967_43_000000000-A3JHG_1_2108_14707_9807")
trash_codes <- c("example", "removing", "sequences")

xdev_remove_sequences(data = data, sequence_names = seqs_to_remove,
                      trash_tags = trash_codes)
#> miseq_sop:
#> 
#>             starts ends   nbases ambigs polymers numns numseqs
#> Minimum:         1  375 249.0000      0 3.000000     0       1
#> 2.5%-tile:       1  375 252.0000      0 4.000000     0    2850
#> 25%-tile:        1  375 252.0000      0 4.000000     0   28491
#> Median:          1  375 253.0000      0 4.000000     0   56981
#> 75%-tile:        1  375 253.0000      0 5.000000     0   85471
#> 97.5%-tile:      1  375 254.0000      0 6.000000     0  111112
#> Maximum:         1  375 256.0000      0 6.000000     0  113961
#> Mean:            1  375 252.7403      0 4.496284     0   56981
#>     Length Overlap_Length Overlap_Start Overlap_End MisMatches Num_Ns
#> 1 250.0000       232.0000      0.000000    248.0000   0.000000      0
#> 2 252.0000       246.0000      1.000000    250.0000   0.000000      0
#> 3 252.0000       249.0000      2.000000    251.0000   0.000000      0
#> 4 253.0000       249.0000      2.000000    251.0000   1.000000      0
#> 5 253.0000       250.0000      2.000000    251.0000   5.000000      0
#> 6 254.0000       251.0000      4.000000    253.0000  26.000000      0
#> 7 270.0000       255.0000     22.000000    256.0000 120.000000      0
#> 8 252.7572       249.1499      2.005367    251.1552   5.161024      0
#>   Expected_Errors
#> 1      1.00000000
#> 2      1.00000000
#> 3      1.00000000
#> 4      1.00000000
#> 5      1.00000000
#> 6      1.00000000
#> 7      4.00000000
#> 8      0.07381921
#>       type trash_code unique total
#> 1 sequence    example      1     1
#> 2 sequence   removing      1     1
#> 3 sequence  sequences      1     1
#> 4      asv    example      1     1
#> 5      asv   removing      1     1
#> 6      asv  sequences      1     1
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  375    249      0        3     0      1.00
#> 2.5%-tile:       1  375    252      0        4     0   2850.00
#> 25%-tile:        1  375    252      0        4     0  28491.00
#> Median:          1  375    253      0        4     0  56981.00
#> 75%-tile:        1  375    253      0        5     0  85471.00
#> 97.5%-tile:      1  375    254      0        6     0 111112.00
#> Maximum:         1  375    256      0        6     0 113961.00
#> Mean:            1  375    252      0        4     0  56981.00
#> scrap_summary:
#>       type trash_code unique total
#> 1 sequence    example      1     1
#> 2 sequence   removing      1     1
#> 3 sequence  sequences      1     1
#> 4      asv    example      1     1
#> 5      asv   removing      1     1
#> 6      asv  sequences      1     1
#> 
#> Number of unique seqs: 2422 
#> Total number of seqs: 113960 
#> 
#> Total number of samples: 19 
#> Total number of treatments: 2 
#> Total number of otus: 531 
#> Total number of otu bin classifications: 531 
#> Total number of asvs: 2422 
#> Total number of asv bin classifications: 2422 
#> Total number of phylotypes: 63 
#> Total number of phylotype bin classifications: 63 
#> Total number of sequence classifications: 2422 
#> Total number of resource references: 2 
#> Total number of custom reports: 1 
#> Your dataset includes metadata 
#> 

# If you look at the scrap report, you the sequences names, listed with the
# trash codes set to "example", "removing", "sequences".

report(data = data, type = "sequence_scrap")
#>                                             id trash_code
#> 1 M00967_43_000000000-A3JHG_1_2101_16474_12783    example
#> 2  M00967_43_000000000-A3JHG_1_1113_12711_3318   removing
#> 3  M00967_43_000000000-A3JHG_1_2108_14707_9807  sequences

# You can see from the get_num_sequences function that the removed
# sequence's abundances are removed from the dataset.

count(data = data, type = "sequence")
#> [1] 113960
```
