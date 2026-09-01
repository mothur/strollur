# Assign samples distances in a [strollur object](https://mothur.org/strollur/reference/strollur.html)

Assign samples distances in a [strollur
object](https://mothur.org/strollur/reference/strollur.html)

## Usage

``` r
xdev_assign_sample_distances(data, table, reference = NULL, verbose = TRUE)
```

## Arguments

- data:

  a [strollur
  object](https://mothur.org/strollur/reference/strollur.html)

- table:

  a 3 column data.frame (sample1 sample2 distance) containing distances
  between your samples

- reference:

  a list created by the function
  [new_reference](https://mothur.org/strollur/reference/new_reference.md).
  Optional.

- verbose:

  a boolean indicating whether or not you want progress messages.
  Default = TRUE.

## Value

an updated [strollur
object](https://mothur.org/strollur/reference/strollur.html)

## Examples

``` r

shared_file <- strollur_example("final.opti_mcc.shared")
dist_file <- strollur_example("final.opti_mcc.jclass.0.03.column.dist")
reference <- strollur::new_reference(name = "jclass estimator distances",
                                     vendor = "mothur_v1.48.6")

data <- strollur::new_dataset("my_dataset")
df <- read_mothur_shared(shared_file)
xdev_assign_bins(data, table = df, bin_type = "otu")
#> Assigned 531 otu bins.
#> my_dataset:
#> 
#> 
#> Number of unique seqs: 531 
#> Total number of seqs: 113963 
#> 
#> Total number of samples: 19 
#> Total number of otus: 531 
#> 

sample_dists <- readr::read_table(dist_file,
                                  col_names = FALSE,
                                  show_col_types = FALSE)
xdev_assign_sample_distances(data, table = sample_dists,
                             reference = reference)
#> Assigned 171 samples distances.
#> Added 1 resource references.
#> my_dataset:
#> 
#> 
#> Number of unique seqs: 531 
#> Total number of seqs: 113963 
#> 
#> Total number of samples: 19 
#> Total number of otus: 531 
#> Total number of resource references: 1 
#> 
data
#> my_dataset:
#> 
#> 
#> Number of unique seqs: 531 
#> Total number of seqs: 113963 
#> 
#> Total number of samples: 19 
#> Total number of otus: 531 
#> Total number of resource references: 1 
#> 
```
