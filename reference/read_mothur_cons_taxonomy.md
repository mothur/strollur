# read_mothur_cons_taxonomy

Read a mothur formatted [cons_taxonomy
file](https://mothur.org/wiki/constaxonomy_file/)

## Usage

``` r
read_mothur_cons_taxonomy(taxonomy)
```

## Arguments

- taxonomy:

  file name, a mothur [consensus taxonomy
  file](https://mothur.org/wiki/constaxonomy_file/). The cons_taxonomy
  file is created by
  [classify.otu](https://mothur.org/wiki/classify.otu/).

## Value

A data.frame containing the bin names, bin abundances and bin
taxonomies.

## Examples

``` r
# You can add the otu assignments and bin taxonomies to the your data set
# using the following:

# read mothur's consensus taxonomy file into a data.frame
otu_data <- read_mothur_cons_taxonomy(strollur_example(
  "final.cons.taxonomy"
))

data <- new_dataset()

# assign abundance only 'otu' bins
assign(data = data, table = otu_data, type = "bins", bin_type = "otu")
#> ℹ Assigned 531 otu bins.
#> [1] 531

# assign consensus taxonomies to 'otu' bins
assign(
  data = data, table = otu_data,
  type = "bin_taxonomy", bin_type = "otu"
)
#> ℹ Assigned 531 otu bin taxonomies.
#> [1] 531
```
