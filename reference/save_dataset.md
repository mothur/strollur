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

## See also

\[load_dataset()\]

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
#> Added 2425 sequences.
#> Assigned 2425 sequence abundances.
#> Assigned 2425 sequence taxonomies.
#> Assigned 531 otu bins.
#> Assigned 19 samples to treatments.

save_dataset(data, "miseq_sop.rds")
#> [1] "miseq_sop.rds"
```
