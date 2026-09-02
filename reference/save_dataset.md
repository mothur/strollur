# Save the [strollur object](https://mothur.org/strollur/reference/strollur.html) to file.

The save_dataset function will save the [strollur
object](https://mothur.org/strollur/reference/strollur.html) to file.

## Usage

``` r
save_dataset(data, file)
```

## Arguments

- data:

  a [strollur
  object](https://mothur.org/strollur/reference/strollur.html)

- file:

  a string containing the file name.

## Value

A file containing the
[`strollur::strollur`](https://mothur.org/strollur/reference/strollur.md)
object

## See also

[`load_dataset()`](https://mothur.org/strollur/reference/load_dataset.md)

## Examples

``` r

data <- strollur::read_mothur(
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

file_name <- file.path(tempdir(), "miseq_sop.rds")
strollur::save_dataset(data, file = file_name)
#> [1] "/tmp/RtmpsPrSxR/miseq_sop.rds"
```
