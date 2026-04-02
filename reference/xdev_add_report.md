# xdev_add_report

Add a report to a
[strollur](https://mothur.org/strollur/reference/strollur.md) object

## Usage

``` r
xdev_add_report(
  data,
  table,
  type = "metadata",
  sequence_name = "sequence_names",
  verbose = TRUE
)
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.md) object

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

## Examples

``` r
# To add a custom report including your contigs assembly data

data <- new_dataset("just for fun", 2)
contigs_report <- readr::read_tsv(strollur_example("final.contigs_report.gz"),
   col_names = TRUE, show_col_types = FALSE)

xdev_add_report(data, contigs_report, "contigs_report", "Name")
#> ℹ Added a contigs_report.

# To add metadata related to your study

metadata <- readr::read_tsv(strollur_example("mouse.dpw.metadata"),
                            col_names = TRUE, show_col_types = FALSE)

xdev_add_report(data, metadata, "metadata")
#> ℹ Added metadata.
```
