# read_mothur_taxonomy

Read a mothur formatted [taxonomy
file](https://mothur.org/wiki/taxonomy_file/)

## Usage

``` r
read_mothur_taxonomy(taxonomy)
```

## Arguments

- taxonomy:

  file name. a mothur [taxonomy
  file](https://mothur.org/wiki/taxonomy_file/), created by
  [classify.seqs](https://mothur.org/wiki/classify.seqs/)

## Value

A data.frame containing the sequences names and sequences taxonomies.

## Examples

``` r
# You can add the sequences and their taxonomies to the your data set
# using the following:

# read mothur's taxonomy file into a data.frame
classification_data <- read_mothur_taxonomy(strollur_example(
  "final.taxonomy.gz"
))

# create a new empty `strollur` object
data <- new_dataset()

# assign sequence classifications
assign(data = data, table = classification_data, type = "sequence_taxonomy")
#> Assigned 2425 sequence taxonomies.
#> [1] 2425
```
