# rdataset_example

rdataset comes bundled with some example files in its \`inst/extdata\`
directory. This function make them easy to access them.

## Usage

``` r
rdataset_example(file = NULL)
```

## Arguments

- file:

  Name of file.

## Details

Get path to rdataset example

## Examples

``` r
rdataset_example()
#> [1] "/home/runner/work/_temp/Library/rdataset/extdata"
rdataset_example("final.fasta")
#> [1] "/home/runner/work/_temp/Library/rdataset/extdata/final.fasta"
```
