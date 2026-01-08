# clone

The clone function will create a deep copy of 'dataset' object. Note:
'dataset' is an R6 class with reference semantics and requires a deep
copy.

## Usage

``` r
clone(data)
```

## Arguments

- data:

  a 'dataset' object

## Value

A 'dataset' object

## Examples

``` r
# For dataset's including sequence data:

data <- read_mothur(
  fasta = rdataset_example("final.fasta"),
  count = rdataset_example("final.count_table"),
  taxonomy = rdataset_example("final.taxonomy"),
  design = rdataset_example("mouse.time.design"),
  otu_list = rdataset_example("final.opti_mcc.list"),
  dataset_name = "miseq_sop"
)
#> ℹ Added 2425 sequences.
#> ℹ Assigned 2425 sequence abundances.
#> ℹ Assigned 2425 sequence taxonomies.
#> ℹ Assigned 531 otu bins.
#> ℹ Assigned 19 samples to treatments.

copy_dataset <- clone(data)
copy_dataset
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
#> 
```
