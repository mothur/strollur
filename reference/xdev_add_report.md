# Add a report to a [strollur](https://mothur.org/strollur/reference/strollur.html) object

Add a report to a
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Usage

``` r
xdev_add_report(
  data,
  table,
  type = "metadata",
  sequence_name = "sequence_name",
  verbose = TRUE
)
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object

- table, :

  a data.frame containing your report.

- type, :

  a string containing the type of report. Options include: "metadata"
  and custom report tags. Default = "metadata".

- sequence_name, :

  a string containing the name of the column in 'table' that contains
  the sequence names. This is used for custom reports, metadata does not
  require a sequence_name column. Default column name is
  'sequence_names'.

- verbose, :

  a boolean whether or not you want progress messages. Default = TRUE.

## Value

an updated
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Examples

``` r

# To add a custom report including your contigs assembly data

data <- new_dataset("just for fun", 2)
contigs_report <- readRDS(strollur_example("miseq_contigs_report.rds"))

xdev_add_report(data, contigs_report, "contigs_report", "Name")
#> Added a contigs_report.
#> just for fun:
#> 
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
#> 
#> Number of unique seqs: 2425 
#> Total number of seqs: 2425 
#> 
#> Total number of custom reports: 1 
#> 

# To add metadata related to your study

metadata <- readRDS(strollur_example("miseq_metadata.rds"))

xdev_add_report(data, metadata, "metadata")
#> Added metadata.
#> just for fun:
#> 
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
#> 
#> Number of unique seqs: 2425 
#> Total number of seqs: 2425 
#> 
#> Total number of custom reports: 1 
#> Your dataset includes metadata 
#> 
```
