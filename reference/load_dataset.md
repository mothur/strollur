# Load strollur object from .rds file

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

## See also

\[save_dataset()\]

## Examples

``` r

data <- load_dataset(strollur_example("miseq_sop.rds"))
data
#> miseq_sop:
#> 
#>             starts ends   nbases ambigs polymers numns    numseqs
#> Minimum:         1  375 249.0000      0 3.000000     0      1.000
#> 2.5%-tile:       1  375 252.0000      0 4.000000     0   2850.075
#> 25%-tile:        1  375 252.0000      0 4.000000     0  28491.750
#> Median:          1  375 253.0000      0 4.000000     0  56982.500
#> 75%-tile:        1  375 253.0000      0 5.000000     0  85473.250
#> 97.5%-tile:      1  375 254.0000      0 6.000000     0 111114.925
#> Maximum:         1  375 256.0000      0 6.000000     0 113964.000
#> Mean:            1  375 252.7406      0 4.496082     0  56982.500
#>     Length Overlap_Length Overlap_Start Overlap_End MisMatches Num_Ns
#> 1 250.0000       232.0000      0.000000    248.0000   0.000000      0
#> 2 252.0000       246.0000      1.000000    250.0000   0.000000      0
#> 3 252.0000       249.0000      2.000000    251.0000   0.000000      0
#> 4 253.0000       249.0000      2.000000    251.0000   1.000000      0
#> 5 253.0000       250.0000      2.000000    251.0000   5.000000      0
#> 6 254.0000       251.0000      4.000000    253.0000  26.000000      0
#> 7 270.0000       255.0000     22.000000    256.0000 120.000000      0
#> 8 252.7575       249.1501      2.005361    251.1555   5.162474      0
#>   Expected_Errors
#> 1      1.00000000
#> 2      1.00000000
#> 3      1.00000000
#> 4      1.00000000
#> 5      1.00000000
#> 6      1.00000000
#> 7      4.00000000
#> 8      0.07385095
#> data frame with 0 columns and 0 rows
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  375    249      0        3     0      1.00
#> 2.5%-tile:       1  375    252      0        4     0   2850.08
#> 25%-tile:        1  375    252      0        4     0  28491.75
#> Median:          1  375    253      0        4     0  56982.50
#> 75%-tile:        1  375    253      0        5     0  85473.25
#> 97.5%-tile:      1  375    254      0        6     0 111114.93
#> Maximum:         1  375    256      0        6     0 113964.00
#> Mean:            1  375    252      0        4     0  56982.50
#> 
#> Number of unique seqs: 2425 
#> Total number of seqs: 113963 
#> 
#> Total number of samples: 19 
#> Total number of treatments: 2 
#> Total number of otus: 531 
#> Total number of otu bin classifications: 531 
#> Total number of asvs: 2425 
#> Total number of asv bin classifications: 2425 
#> Total number of phylotypes: 63 
#> Total number of phylotype bin classifications: 63 
#> Total number of sequence classifications: 2425 
#> Total number of resource references: 2 
#> Total number of custom reports: 1 
#> Your dataset includes metadata 
#> 
```
