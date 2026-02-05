# xdev_summarize

Summarize the sequences data, custom reports, and scrapped data in a
[dataset](https://mothur.org/strollur/reference/dataset.md) object

## Usage

``` r
xdev_summarize(data, type = "sequences", report_type = NULL)
```

## Arguments

- data, :

  a [dataset](https://mothur.org/strollur/reference/dataset.md) object

- type, :

  string containing the type of data you want the number of. Options
  include: "sequences", "reports" and "scrap". Default = "sequences".

- report_type, :

  string containing the report type you would summarized. For example,
  the miseq_sop_example includes contigs assembly data and can be
  accessed with report_type = "contigs_report". Default = NULL.

## Value

data.frame()

## Examples

``` r
 data <- miseq_sop_example()
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

 # summarize FASTA data
 xdev_summarize(data = data, type = "sequences")
#>             starts ends   nbases ambigs polymers numns    numseqs
#> Minimum:         1  375 249.0000      0 3.000000     0      1.000
#> 2.5%-tile:       1  375 252.0000      0 3.000000     0   2850.075
#> 25%-tile:        1  375 252.0000      0 4.000000     0  28491.750
#> Median:          1  375 252.0000      0 4.000000     0  56982.500
#> 75%-tile:        1  375 253.0000      0 5.000000     0  85473.250
#> 97.5%-tile:      1  375 253.0000      0 6.000000     0 111114.925
#> Maximum:         1  375 256.0000      0 6.000000     0 113963.000
#> Mean:            1  375 252.4472      0 4.368699     0      0.000

 # summarize contigs_report
 xdev_summarize(data = data, type = "reports",
                 report_type = "contigs_report")
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

 # remove sample 'F3D0'
 xdev_remove_samples(data = data, samples = c("F3D0"))

 # summarize FASTA data after removal of sample F3D0
 xdev_summarize(data = data, type = "sequences")
#>             starts ends   nbases ambigs polymers numns  numseqs
#> Minimum:         1  375 249.0000      0 3.000000     0      1.0
#> 2.5%-tile:       1  375 252.0000      0 3.000000     0   2695.3
#> 25%-tile:        1  375 252.0000      0 4.000000     0  26944.0
#> Median:          1  375 252.0000      0 4.000000     0  53887.0
#> 75%-tile:        1  375 253.0000      0 5.000000     0  80830.0
#> 97.5%-tile:      1  375 253.0000      0 6.000000     0 105078.7
#> Maximum:         1  375 256.0000      0 6.000000     0 107772.0
#> Mean:            1  375 252.4444      0 4.361532     0      0.0

 # summarize scrapped data
 xdev_summarize(data = data, type = "scrap")
#>        type     trash_code unique total
#> 1  sequence remove_samples    101   109
#> 2       otu remove_samples     14    14
#> 3       asv remove_samples    101   109
#> 4 phylotype remove_samples      2     2
```
