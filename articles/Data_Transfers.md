# Data_Transfers

The *rdataset* package stores the data associated with your microbial
DNA analysis. This tutorial will explain how to save, load, copy,
export, and import your dataset object. If you haven’t reviewed the
[Getting Started](Getting-Started.md) tuturial, we recommend you start
there. First let’s load the rdataset package.

``` r
library(rdataset)
#> 
#> Attaching package: 'rdataset'
#> The following objects are masked from 'package:base':
#> 
#>     assign, names, summary
```

We can use the **miseq_sop_example** function to create a dataset object
from the [Miseq SOP Example](https://mothur.org/wiki/miseq_sop/).

``` r
miseq <- miseq_sop_example()
#> ℹ Added 2425 sequences.
#> ℹ Assigned 2425 sequence abundances.
#> ℹ Assigned 2425 sequence taxonomies.
#> ℹ Assigned 531 otu bins.
#> ℹ Assigned 2425 asv bins.
#> ℹ Assigned 63 phylotype bins.
#> ℹ Assigned 19 samples to treatments.
#> ℹ Assigned 531 otu bin taxonomies.
#> ℹ Assigned 531 otu bin representative sequences.
#> ℹ Added metadata.
#> ℹ Added 2 resource references.
#> ℹ Added a contigs_report.
miseq
#> miseq_sop:
#> 
#> sequence_summary:
#>             starts ends   nbases ambigs polymers numns    numseqs
#> Minimum:         1  375 249.0000      0 3.000000     0      1.000
#> 2.5%-tile:       1  375 252.0000      0 3.000000     0   2850.075
#> 25%-tile:        1  375 252.0000      0 4.000000     0  28491.750
#> Median:          1  375 252.0000      0 4.000000     0  56982.500
#> 75%-tile:        1  375 253.0000      0 5.000000     0  85473.250
#> 97.5%-tile:      1  375 253.0000      0 6.000000     0 111114.925
#> Maximum:         1  375 256.0000      0 6.000000     0 113963.000
#> Mean:            1  375 252.4472      0 4.368699     0      0.000
#> Unique seqs:  2425 
#> Total seqs:   113963 
#> 
#> contigs_report :
#>             Expected_Errors   Length MisMatches Num_Ns Overlap_End
#> Minimum:       0.0000452496 250.0000   0.000000      0    248.0000
#> 2.5%-tile:     0.0016101600 252.0000   0.000000      0    251.0000
#> 25%-tile:      0.0028177700 252.0000   0.000000      0    251.0000
#> Median:        0.0062948698 252.0000   2.000000      0    251.0000
#> 75%-tile:      0.0264780000 253.0000   4.000000      0    251.0000
#> 97.5%-tile:    1.1412299871 253.0000  76.000000      0    252.0000
#> Maximum:       3.0126200000 270.0000 120.000000      0    256.0000
#> Mean:          0.0984788569 252.5128   7.534147      0    251.0762
#>             Overlap_Length Overlap_Start
#> Minimum:          232.0000      0.000000
#> 2.5%-tile:        248.0000      1.000000
#> 25%-tile:         249.0000      2.000000
#> Median:           249.0000      2.000000
#> 75%-tile:         249.0000      2.000000
#> 97.5%-tile:       251.0000      3.000000
#> Maximum:          255.0000     22.000000
#> Mean:             249.2136      1.862692
#> Unique seqs:  2425 
#> Total seqs:   113963 
#> 
#> 
#> Sample   Total:
#> F3D0 6191 
#> F3D1 4652 
#> F3D141   4656 
#> F3D142   2423 
#> F3D143   2403 
#> F3D144   3449 
#> F3D145   5532 
#> F3D146   3831 
#> F3D147   12430 
#> F3D148   9465 
#> F3D149   10014 
#> F3D150   4126 
#> F3D2 15686 
#> F3D3 5199 
#> F3D5 3469 
#> F3D6 6394 
#> F3D7 4055 
#> F3D8 4253 
#> F3D9 5735 
#> 
#> Treatment   Total:
#> Early    55634 
#> Late 58329 
#> 
#> Number of unique seqs: 2425 
#> Total number of seqs: 113963 
#> Total number of otus: 531 
#> Total number of asvs: 2425 
#> Total number of phylotypes: 63
```

## Save

## Load

## Copy

## Export

## Import
