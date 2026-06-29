# xdev_assign_bin_representative_sequences

Assign representative sequences to bins.

## Usage

``` r
xdev_assign_bin_representative_sequences(
  data,
  table,
  bin_type = "otu",
  reference = NULL,
  bin_name = "bin_name",
  sequence_name = "sequence_name",
  verbose = TRUE
)
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object

- table, :

  a data.frame containing bin representative assignments

- bin_type:

  a string indicating the type of bin assignments. Default "otu".

- reference, :

  a list created by the function \[new_reference\]. Optional.

- bin_name, :

  a string containing the name of the column in 'table' that contains
  the bin names. Default column name is 'bin_name'.

- sequence_name:

  a string containing the name of the column in 'table' that contains
  the bin names. Default column name is 'sequence_name'.

- verbose, :

  a boolean whether or not you want progress messages. Default = TRUE.

## Value

an updated
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Examples

``` r

  miseq <- miseq_sop_example()
#> Added 2425 sequences.
#> Assigned 2425 sequence abundances.
#> Assigned 2425 sequence taxonomies.
#> Assigned 531 otu bins.
#> Assigned 2425 asv bins.
#> Assigned 63 phylotype bins.
#> Assigned 19 samples to treatments.
#> Assigned 531 otu bin taxonomies.
#> Assigned 531 otu bin representative sequences.
#> Added a metadata report.
#> Added 2 resource references.
#> Added a contigs_report report.

  bin_reps <- readRDS(strollur_example(
         "miseq_representative_sequences.rds"))

  xdev_assign_bin_representative_sequences(data = miseq,
                                      table = bin_reps,
                                      bin_type = "otu")
#> Assigned 531 otu bin representative sequences.
#> miseq_sop:
#> 
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  375    249      0        3     0      1.00
#> 2.5%-tile:       1  375    252      0        4     0   2849.08
#> 25%-tile:        1  375    252      0        4     0  28490.75
#> Median:          1  375    253      0        4     0  56981.50
#> 75%-tile:        1  375    253      0        5     0  85472.25
#> 97.5%-tile:      1  375    254      0        6     0 111113.93
#> Maximum:         1  375    256      0        6     0 113963.00
#> Mean:            1  375    252      0        4     0  56981.64
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
#> Total number of custom reports: 2 
#> 
```
