# strollur_example

strollur comes bundled with some example files in its \`inst/extdata\`
directory. This function make them easy to access them.

## Usage

``` r
strollur_example(file = NULL)
```

## Arguments

- file:

  Name of file.

## Details

Get path to strollur example

## Examples

``` r
strollur_example()
#> [1] "/home/runner/work/_temp/Library/strollur/extdata"
strollur_example("final.fasta.gz")
#> [1] "/home/runner/work/_temp/Library/strollur/extdata/final.fasta.gz"
```
