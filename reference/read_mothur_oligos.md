# read_mothur_oligos

Read a mothur formatted [oligos
file](https://mothur.org/wiki/oligos_file/)

## Usage

``` r
read_mothur_oligos(oligos)
```

## Arguments

- oligos:

  file name. a mothur formatted [oligos
  file](https://mothur.org/wiki/oligos_file/)

## Value

A data.frame containing the oligos data.

## Examples

``` r

oligos <- read_mothur_oligos(strollur_example("paired_read.oligos"))

# Create a new dataset and add your oligos data

data <- new_dataset() |>
  add(
    table = oligos,
    type = "report",
    report_type = "paired_oligos"
  )
#> Added a paired_oligos report.
```
