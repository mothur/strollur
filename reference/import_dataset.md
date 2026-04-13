# Import strollur object from exported data.frame.

The import_dataset function will create a
[strollur](https://mothur.org/strollur/reference/strollur.html) object
from the exported table of a
[strollur](https://mothur.org/strollur/reference/strollur.html) object.

## Usage

``` r
import_dataset(table)
```

## Arguments

- table:

  a table containing the data from a
  [strollur](https://mothur.org/strollur/reference/strollur.html)
  object. You can create the table using 'export(data)'.

## Value

a [strollur](https://mothur.org/strollur/reference/strollur.html) object

## See also

\[export_dataset()\]

## Examples

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
data <- import_dataset(export_dataset(miseq))
#> ℹ Added 2425 sequences.
#> ℹ Assigned 2425 sequence taxonomies.
#> ℹ Assigned 2425 sequence abundances.
#> ℹ Assigned 531 otu bins.
#> ℹ Assigned 531 otu bin representative sequences.
#> ℹ Assigned 531 otu bin taxonomies.
#> ℹ Assigned 2425 asv bins.
#> ℹ Assigned 2425 asv bin taxonomies.
#> ℹ Assigned 63 phylotype bins.
#> ℹ Assigned 63 phylotype bin taxonomies.
#> ℹ Added metadata.
#> ℹ Added 2 resource references.
#> ℹ Added a contigs_report.
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
