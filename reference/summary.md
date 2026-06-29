# Summarize the sequences data, custom reports, and scrapped data in a [strollur](https://mothur.org/strollur/reference/strollur.html) object

Summarize the sequences data, custom reports, and scrapped data in a
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Usage

``` r
summary(data, type = "sequence", report_type = NULL, verbose = TRUE)
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object

- type, :

  string containing the type of data you want the number of. Options
  include: "sequence", "report" and "scrap". Default = "sequence".

- report_type, :

  string containing the report type you would summarized. For example,
  the miseq_sop_example includes contigs assembly data and can be
  accessed with report_type = "contigs_report". Default = NULL.

- verbose, :

  boolean indicating whether or not you want progress messages. Default
  = TRUE.

## Value

data.frame

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

# To get the summary of your FASTA data
summary(data = miseq, type = "sequence")
#>             starts ends   nbases ambigs polymers numns    numseqs
#> Minimum:         1  375 249.0000      0 3.000000     0      1.000
#> 2.5%-tile:       1  375 252.0000      0 4.000000     0   2849.075
#> 25%-tile:        1  375 252.0000      0 4.000000     0  28490.750
#> Median:          1  375 253.0000      0 4.000000     0  56981.500
#> 75%-tile:        1  375 253.0000      0 5.000000     0  85472.250
#> 97.5%-tile:      1  375 254.0000      0 6.000000     0 111113.925
#> Maximum:         1  375 256.0000      0 6.000000     0 113963.000
#> Mean:            1  375 252.7406      0 4.496082     0  56981.643
#>             starts ends   nbases ambigs polymers numns    numseqs
#> Minimum:         1  375 249.0000      0 3.000000     0      1.000
#> 2.5%-tile:       1  375 252.0000      0 4.000000     0   2849.075
#> 25%-tile:        1  375 252.0000      0 4.000000     0  28490.750
#> Median:          1  375 253.0000      0 4.000000     0  56981.500
#> 75%-tile:        1  375 253.0000      0 5.000000     0  85472.250
#> 97.5%-tile:      1  375 254.0000      0 6.000000     0 111113.925
#> Maximum:         1  375 256.0000      0 6.000000     0 113963.000
#> Mean:            1  375 252.7406      0 4.496082     0  56981.643

# summarize contigs_report
summary(data = miseq, type = "report", report_type = "contigs_report")
#>               Length Overlap_Length Overlap_Start Overlap_End MisMatches Num_Ns
#> Minimum:    250.0000       232.0000      0.000000    248.0000   0.000000      0
#> 2.5%-tile:  252.0000       246.0000      1.000000    250.0000   0.000000      0
#> 25%-tile:   252.0000       249.0000      2.000000    251.0000   0.000000      0
#> Median:     253.0000       249.0000      2.000000    251.0000   1.000000      0
#> 75%-tile:   253.0000       250.0000      2.000000    251.0000   5.000000      0
#> 97.5%-tile: 254.0000       251.0000      4.000000    253.0000  26.000000      0
#> Maximum:    270.0000       255.0000     22.000000    256.0000 120.000000      0
#> Mean:       252.7575       249.1501      2.005361    251.1555   5.162474      0
#>             Expected_Errors
#> Minimum:         1.00000000
#> 2.5%-tile:       1.00000000
#> 25%-tile:        1.00000000
#> Median:          1.00000000
#> 75%-tile:        1.00000000
#> 97.5%-tile:      1.00000000
#> Maximum:         4.00000000
#> Mean:            0.07385095
#>               Length Overlap_Length Overlap_Start Overlap_End MisMatches Num_Ns
#> Minimum:    250.0000       232.0000      0.000000    248.0000   0.000000      0
#> 2.5%-tile:  252.0000       246.0000      1.000000    250.0000   0.000000      0
#> 25%-tile:   252.0000       249.0000      2.000000    251.0000   0.000000      0
#> Median:     253.0000       249.0000      2.000000    251.0000   1.000000      0
#> 75%-tile:   253.0000       250.0000      2.000000    251.0000   5.000000      0
#> 97.5%-tile: 254.0000       251.0000      4.000000    253.0000  26.000000      0
#> Maximum:    270.0000       255.0000     22.000000    256.0000 120.000000      0
#> Mean:       252.7575       249.1501      2.005361    251.1555   5.162474      0
#>             Expected_Errors
#> Minimum:         1.00000000
#> 2.5%-tile:       1.00000000
#> 25%-tile:        1.00000000
#> Median:          1.00000000
#> 75%-tile:        1.00000000
#> 97.5%-tile:      1.00000000
#> Maximum:         4.00000000
#> Mean:            0.07385095

# remove sample 'F3D0' to produce a scrap report
xdev_remove_samples(data = miseq, samples = c("F3D0"))
#> miseq_sop:
#> 
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  375    249      0        3     0      1.00
#> 2.5%-tile:       1  375    252      0        4     0   2694.30
#> 25%-tile:        1  375    252      0        4     0  26943.00
#> Median:          1  375    253      0        4     0  53886.00
#> 75%-tile:        1  375    253      0        5     0  80829.00
#> 97.5%-tile:      1  375    254      0        6     0 105077.70
#> Maximum:         1  375    256      0        6     0 107772.00
#> Mean:            1  375    252      0        4     0  53886.14
#> 
#> scrap_summary:
#>        type     trash_code unique total
#> 1  sequence remove_samples    101   109
#> 2       otu remove_samples     14    14
#> 3       asv remove_samples    101   109
#> 4 phylotype remove_samples      2     2
#> 
#> Number of unique seqs: 2324 
#> Total number of seqs: 107772 
#> 
#> Total number of samples: 18 
#> Total number of treatments: 2 
#> Total number of otus: 517 
#> Total number of otu bin classifications: 517 
#> Total number of asvs: 2324 
#> Total number of asv bin classifications: 2324 
#> Total number of phylotypes: 61 
#> Total number of phylotype bin classifications: 61 
#> Total number of sequence classifications: 2324 
#> Total number of resource references: 2 
#> Total number of custom reports: 2 
#> 

# summarize FASTA data after removal of sample F3D0
summary(data = miseq, type = "sequence")
#>             starts ends   nbases ambigs polymers numns   numseqs
#> Minimum:         1  375 249.0000      0 3.000000     0      1.00
#> 2.5%-tile:       1  375 252.0000      0 4.000000     0   2694.30
#> 25%-tile:        1  375 252.0000      0 4.000000     0  26943.00
#> Median:          1  375 253.0000      0 4.000000     0  53886.00
#> 75%-tile:        1  375 253.0000      0 5.000000     0  80829.00
#> 97.5%-tile:      1  375 254.0000      0 6.000000     0 105077.70
#> Maximum:         1  375 256.0000      0 6.000000     0 107772.00
#> Mean:            1  375 252.7345      0 4.493546     0  53886.14
#>             starts ends   nbases ambigs polymers numns   numseqs
#> Minimum:         1  375 249.0000      0 3.000000     0      1.00
#> 2.5%-tile:       1  375 252.0000      0 4.000000     0   2694.30
#> 25%-tile:        1  375 252.0000      0 4.000000     0  26943.00
#> Median:          1  375 253.0000      0 4.000000     0  53886.00
#> 75%-tile:        1  375 253.0000      0 5.000000     0  80829.00
#> 97.5%-tile:      1  375 254.0000      0 6.000000     0 105077.70
#> Maximum:         1  375 256.0000      0 6.000000     0 107772.00
#> Mean:            1  375 252.7345      0 4.493546     0  53886.14

# summarize scrapped data -
# sequences and bins scrapped by removing the sample "F3D0"
summary(data = miseq, type = "scrap")
#>        type     trash_code unique total
#> 1  sequence remove_samples    101   109
#> 2       otu remove_samples     14    14
#> 3       asv remove_samples    101   109
#> 4 phylotype remove_samples      2     2
#>        type     trash_code unique total
#> 1  sequence remove_samples    101   109
#> 2       otu remove_samples     14    14
#> 3       asv remove_samples    101   109
#> 4 phylotype remove_samples      2     2
```
