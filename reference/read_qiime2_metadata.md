# read_qiime2_metadata

Read a [qiime2](https://qiime2.org) .tsv table containing metadata.

## Usage

``` r
read_qiime2_metadata(metadata)
```

## Arguments

- metadata:

  file name, a qiime2 .tsv file containing metadata about your analysis.

## Value

A data.frame containing metadata

## Examples

``` r
metadata <- read_qiime2_metadata(strollur_example(
  "sample_metadata.tsv"
))
```
