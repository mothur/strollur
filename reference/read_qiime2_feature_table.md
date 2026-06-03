# read_qiime2_feature_table

Read a [qiime2](https://qiime2.org) qza containing bin data

## Usage

``` r
read_qiime2_feature_table(
  qza,
  dir_path = NULL,
  remove_unpacked_artifacts = TRUE
)
```

## Arguments

- qza:

  file name, a qiime2 .qza file containing bin data.

- dir_path:

  a string containing the name of directory where the artifacts files
  should be unpacked. Default = current working directory.

- remove_unpacked_artifacts:

  boolean, When TRUE, the artifact's temporary directories will be
  removed after processing. Default = TRUE.

## Value

A list containing artifact

## Examples

``` r

if (requireNamespace("h5lite", quietly = TRUE)) {
  artifact <- read_qiime2_feature_table(strollur_example("table.qza"))

  # access the bin assignment table

  artifact$data

  # to create a `strollur` object with your data

  data <- new_dataset("my_data")

  assign(data = data, table = artifact$data, type = "bin")
  data
} else {
  message(paste(
    "To use this functionality you have to install the",
    "h5lite package."
  ))
}
#> Assigned 759 otu bins.
#> my_data:
#> 
#> 
#> Number of unique seqs: 759 
#> Total number of seqs: 157298 
#> 
#> Total number of samples: 34 
#> Total number of otus: 759 
#> 
```
