# xdev_remove_lineages

Designed with package integration in mind, the remove lineages function
allows you to remove contaminents from a
[strollur](https://mothur.org/strollur/reference/strollur.html)

## Usage

``` r
xdev_remove_lineages(data, contaminants, reason = "contaminant")
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object.

- contaminants, :

  vector of strings containing the taxonomies you would like to remove

- reason, :

  a string containing reason you are removing the lineages. Default =
  "contaminant".

## Value

an updated
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Examples

``` r
data <- read_mothur(fasta = strollur_example("final.fasta.gz"),
                       count = strollur_example("final.count_table.gz"),
                       taxonomy = strollur_example("final.taxonomy.gz"),
                       design = strollur_example("mouse.time.design"),
                       otu_list = strollur_example("final.opti_mcc.list.gz"),
                       dataset_name = "miseq_sop")
#> Added 2425 sequences.
#> Assigned 2425 sequence abundances.
#> Assigned 2425 sequence taxonomies.
#> Assigned 531 otu bins.
#> Assigned 19 samples to treatments.

contaminants <- c("Chloroplast", "Mitochondria", "unknown", "Archaea",
 "Eukaryota")

xdev_remove_lineages(data = data, contaminants = contaminants)
#> miseq_sop:
#> 
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  375    249      0        3     0      1.00
#> 2.5%-tile:       1  375    252      0        4     0   2850.08
#> 25%-tile:        1  375    252      0        4     0  28491.75
#> Median:          1  375    253      0        4     0  56982.50
#> 75%-tile:        1  375    253      0        5     0  85473.25
#> 97.5%-tile:      1  375    254      0        6     0 111114.93
#> Maximum:         1  375    256      0        6     0 113964.00
#> Mean:            1  375    252      0        4     0  56982.50
#> 
#> Number of unique seqs: 2425 
#> Total number of seqs: 113963 
#> 
#> Total number of samples: 19 
#> Total number of treatments: 2 
#> Total number of otus: 531 
#> Total number of otu bin classifications: 531 
#> Total number of sequence classifications: 2425 
#> 
```
