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
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  375    249      0        3     0      1.00
#> 2.5%-tile:       1  375    252      0        4     0   2849.00
#> 25%-tile:        1  375    252      0        4     0  28490.00
#> Median:          1  375    253      0        4     0  56980.00
#> 75%-tile:        1  375    253      0        5     0  85470.00
#> 97.5%-tile:      1  375    254      0        6     0 111111.00
#> Maximum:         1  375    256      0        6     0 113960.00
#> Mean:            1  375    252      0        4     0  56980.14
#> 
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
#>                                  sequence_name trash_code
#> 1 M00967_43_000000000-A3JHG_1_2101_16474_12783    example
#> 2  M00967_43_000000000-A3JHG_1_1113_12711_3318   removing
#> 3  M00967_43_000000000-A3JHG_1_2108_14707_9807  sequences

# You can see from the get_num_sequences function that the removed
# sequence's abundances are removed from the dataset.

count(data = data, type = "sequence")
#> [1] 113960
```
