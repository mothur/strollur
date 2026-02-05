# xdev_remove_lineages

Designed with package integration in mind, the remove lineages function
allows you to remove contaminents from a
[dataset](https://mothur.org/strollur/reference/dataset.md)

## Usage

``` r
xdev_remove_lineages(data, contaminants, reason = "contaminant")
```

## Arguments

- data, :

  a [dataset](https://mothur.org/strollur/reference/dataset.md) object.

- contaminants, :

  vector of strings containing the taxonomies you would like to remove

- reason, :

  a string containing reason you are removing the lineages. Default =
  "contaminant".

## Examples

``` r
data <- read_mothur(fasta = strollur_example("final.fasta.gz"),
                       count = strollur_example("final.count_table.gz"),
                       taxonomy = strollur_example("final.taxonomy.gz"),
                       design = strollur_example("mouse.time.design"),
                       otu_list = strollur_example("final.opti_mcc.list.gz"),
                       dataset_name = "miseq_sop")
#> ℹ Added 2425 sequences.
#> ℹ Assigned 2425 sequence abundances.
#> ℹ Assigned 2425 sequence taxonomies.
#> ℹ Assigned 531 otu bins.
#> ℹ Assigned 19 samples to treatments.

contaminants <- c("Chloroplast", "Mitochondria", "unknown", "Archaea",
 "Eukaryota")

xdev_remove_lineages(data = data, contaminants = contaminants)
```
