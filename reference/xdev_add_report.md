# Add a report to a [strollur](https://mothur.org/strollur/reference/strollur.html) object

Add a report to a
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Usage

``` r
xdev_add_report(
  data,
  table,
  type = "report",
  sequence_name = "none",
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

  a string containing the type of report. Default = "report".

- sequence_name, :

  a string. If your report relates to the sequence data,
  \`sequence_name\` should contain the name of the column in 'table'
  that contains the sequence names. Default = 'none'.

- verbose, :

  a boolean whether or not you want progress messages. Default = TRUE.

## Value

an updated
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Examples

``` r

# To add a custom report including your contigs assembly data

data <- new_dataset("just for fun")
contigs_report <- readRDS(strollur_example("miseq_contigs_report.rds"))

xdev_add_report(data, contigs_report, "contigs_report", "Name")
#> Added a contigs_report report.
#> just for fun:
#> 
#> 
#> Number of unique seqs: 2425 
#> Total number of seqs: 2425 
#> 
#> Total number of custom reports: 1 
#> 

# To add metadata related to your study

metadata <- readRDS(strollur_example("miseq_metadata.rds"))

xdev_add_report(data, metadata, "metadata")
#> Added a metadata report.
#> just for fun:
#> 
#> 
#> Number of unique seqs: 2425 
#> Total number of seqs: 2425 
#> 
#> Total number of custom reports: 2 
#> 
```
