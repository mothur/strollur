# Add a report to a [strollur object](https://mothur.org/strollur/reference/strollur.html)

Add a report to a [strollur
object](https://mothur.org/strollur/reference/strollur.html)

## Usage

``` r
xdev_add_report(
  data,
  table,
  reference = NULL,
  type = "report",
  sequence_name = "none",
  verbose = TRUE
)
```

## Arguments

- data:

  a [strollur
  object](https://mothur.org/strollur/reference/strollur.html)

- table:

  a data.frame containing your report.

- reference:

  a list created by the function
  [new_reference](https://mothur.org/strollur/reference/new_reference.md).
  Optional.

- type:

  a string containing the type of report. Default = "report".

- sequence_name, :

  a string. If your report relates to the sequence data, `sequence_name`
  should contain the name of the column in 'table' that contains the
  sequence names. Default = 'none'.

- verbose:

  a boolean whether or not you want progress messages. Default = TRUE.

## Value

an updated [strollur
object](https://mothur.org/strollur/reference/strollur.html)

## Examples

``` r

# To add a custom report including your contigs assembly data

data <- strollur::new_dataset("just for fun")
contigs_report <- readRDS(strollur_example("miseq_contigs_report.rds"))

strollur::xdev_add_report(data,
                            table = contigs_report,
                            type = "contigs_report",
                            sequence_name = "Name")
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

strollur::xdev_add_report(data,
                            table = metadata,
                            type = "metadata")
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
