# save_dataset

The save_dataset function will save the [dataset](dataset.md) object to
file.

## Usage

``` r
save_dataset(dataset, file)
```

## Arguments

- dataset:

  a [dataset](dataset.md) object

- file:

  a string containing the file name.

## Value

A file containing the 'dataset' object

## Examples

``` r
dataset <- read_mothur(
  fasta = rdataset_example("final.fasta"),
  count = rdataset_example("final.count_table"),
  taxonomy = rdataset_example("final.taxonomy"),
  design = rdataset_example("mouse.time.design"),
  otu_list = rdataset_example("final.opti_mcc.list"),
  dataset_name = "miseq_sop"
)
#> ℹ Added 2425 sequences.
#> ℹ Assigned 2425 sequence abundances.
#> ℹ Assigned 2425 sequence taxonomies.
#> ℹ Assigned 531 otu bins.
#> ℹ Assigned 19 samples to treatments.

save_dataset(dataset, "miseq_sop.rds")
#> [1] "miseq_sop.rds"
```
