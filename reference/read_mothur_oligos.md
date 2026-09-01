# Read a mothur formatted [oligos file](https://mothur.org/wiki/oligos_file/)

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

oligos <-
     strollur::read_mothur_oligos(strollur_example("paired_read.oligos"))

# Create a new dataset and add your oligos data

data <- strollur::new_dataset() |>
  strollur::add(
    table = oligos,
    type = "report",
    report_type = "oligos"
  )
#> Added an oligos report.
```
