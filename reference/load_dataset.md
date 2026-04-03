# load_dataset

The load_dataset function will create a
[strollur](https://mothur.org/strollur/reference/strollur.html) object
from an RDS file.

## Usage

``` r
load_dataset(file)
```

## Arguments

- file:

  a string containing the .rds file name.

## Value

a [strollur](https://mothur.org/strollur/reference/strollur.html) object

## Examples

``` r
data <- load_dataset(strollur_example("miseq_sop.rds"))
data
#> miseq_sop:
#> 
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  375    249      0        3     0      1.00
#> 2.5%-tile:       1  375    252      0        3     0   2850.08
#> 25%-tile:        1  375    252      0        4     0  28491.75
#> Median:          1  375    252      0        4     0  56982.50
#> 75%-tile:        1  375    253      0        5     0  85473.25
#> 97.5%-tile:      1  375    253      0        6     0 111114.93
#> Maximum:         1  375    256      0        6     0 113963.00
#> Mean:            1  375    252      0        4     0      0.00
#> 
#> Number of unique seqs: 2425 
#> Total number of seqs: 113963 
#> 
#> Total number of samples: 19 
#> Total number of treatments: 2 
#> Total number of otus: 531 
#> Total number of asvs: 2425 
#> Total number of phylotypes: 63 
#> Total number of resource references: 2 
#> Total number of custom reports: 1 
#> Your dataset includes metadata 
#> 
```
