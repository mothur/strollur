# Create a [strollur object](https://mothur.org/strollur/reference/strollur.html) from a biom formatted file.

The read_biom function reads a [biom formatted](https://biom-format.org)
file and creates a
[`strollur::strollur`](https://mothur.org/strollur/reference/strollur.md)
object.

## Usage

``` r
read_biom(biom)
```

## Arguments

- biom:

  string containing the name of the biom file.

## Value

A [strollur object](https://mothur.org/strollur/reference/strollur.html)

## Examples

``` r

if (requireNamespace("h5lite", quietly = TRUE)) {
  data <- strollur::read_biom(biom = strollur_example("miseq_sop.otu.biom"))
  data
} else {
  message(paste(
    "To use this functionality you have to install the",
    "h5lite package."
  ))
}
#> Added 531 sequences.
#> Assigned 531 sequence abundances.
#> Assigned 531 asv bins.
#> Assigned 531 asv bin taxonomies.
#> Added a report report.
#> miseq_sop:
#> 
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  375    250      0        3     0      1.00
#> 2.5%-tile:       1  375    252      0        4     0   2850.05
#> 25%-tile:        1  375    252      0        4     0  28491.50
#> Median:          1  375    253      0        4     0  56982.00
#> 75%-tile:        1  375    253      0        5     0  85472.50
#> 97.5%-tile:      1  375    254      0        6     0 111113.95
#> Maximum:         1  375    256      0        6     0 113963.00
#> Mean:            1  375    252      0        4     0  56982.00
#> 
#> Number of unique seqs: 531 
#> Total number of seqs: 113963 
#> 
#> Total number of samples: 19 
#> Total number of asvs: 531 
#> Total number of asv bin classifications: 531 
#> Total number of custom reports: 1 
#> 
```
