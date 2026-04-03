# save_dataset

The save_dataset function will save the
[strollur](https://mothur.org/strollur/reference/strollur.html) object
to file.

## Usage

``` r
save_dataset(data, file)
```

## Arguments

- data:

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object

- file:

  a string containing the file name.

## Value

A file containing the \`strollur\` object

## Examples

``` r
data <- read_mothur(
  fasta = strollur_example("final.fasta.gz"),
  count = strollur_example("final.count_table.gz"),
  taxonomy = strollur_example("final.taxonomy.gz"),
  design = strollur_example("mouse.time.design"),
  otu_list = strollur_example("final.opti_mcc.list.gz"),
  dataset_name = "miseq_sop"
)
#> ℹ Added 2425 sequences.
#> ℹ Assigned 2425 sequence abundances.
#> ℹ Assigned 2425 sequence taxonomies.
#> ℹ Assigned 531 otu bins.
#> ℹ Assigned 19 samples to treatments.

save_dataset(data, "miseq_sop.rds")
#> [1] "miseq_sop.rds"
```
