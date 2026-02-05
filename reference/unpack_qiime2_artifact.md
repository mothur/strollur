# unpack_qiime2_artifact

The unpack_qiime2_artifact function reads .qza files created by
[qiime2](https://qiime2.org), and returns the artifact.

To generate the various input files you can follow [qiime
moving-pictures](https://amplicon-docs.qiime2.org/en/latest/tutorials/moving-pictures.html).

## Usage

``` r
unpack_qiime2_artifact(qza, dir_path = NULL)
```

## Arguments

- qza:

  filename, a .qza file containing artifact

- dir_path:

  a string containing the name of directory where the artifacts files
  should be written. Default = current working directory.

## Value

A unpacked qza artifact

## Examples

``` r
# Using the example files from moving-pictures

artifact <- unpack_qiime2_artifact(
  qza = strollur_example("table.qza"),
)
```
