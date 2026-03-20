# add

Add sequences, reports, metadata or resource references to a
[dataset](https://mothur.org/strollur/reference/dataset.md) object

## Usage

``` r
add(
  data,
  table,
  type = "sequences",
  report_type = NULL,
  table_names = list(sequence_name = "sequence_names", sequence = "sequences", comment =
    "comments", reference_name = "reference_names", reference_version =
    "reference_versions", reference_usage = "reference_usages", reference_note =
    "reference_notes", reference_url = "reference_urls"),
  reference = NULL,
  verbose = TRUE
)
```

## Arguments

- data, :

  a [dataset](https://mothur.org/strollur/reference/dataset.md) object

- table, :

  a data.frame containing the data you wish to add.

- type, :

  a string containing the type of data. Options include: 'sequences',
  'references' 'metadata' and 'reports'.

- report_type, :

  a string containing the type of report you are adding. Options
  include: 'metadata' and custom reports.

- table_names, :

  named list used to indicate the names of the columns in the table. By
  default:

  table_names \<- list(sequence_name = "sequence_names", comment =
  "comments", sequence = "sequences", reference_name =
  "reference_names", reference_version = "reference_versions",
  reference_usage = "reference_usages", reference_note =
  "reference_notes", reference_url = "reference_urls")

  In table_names, 'sequence_name' is a string containing the name of the
  column in 'table' that contains the sequence names. It is used when
  you are adding FASTA data. Default column name is 'sequence_names'.

  In table_names, 'sequence' is a string containing the name of the
  column in 'table' that contains the sequence nucleotide strings. It is
  used when you are adding FASTA data. Default column name is
  'sequences'.

  In table_names, 'comment' is a string containing the name of the
  column in 'table' that contains the sequence comments. It is used when
  you are adding FASTA data. Default column name is 'comments'.

  In table_names, 'reference_name' is a string containing the name of
  the column in 'table' that contains the reference names. It is used
  when you are adding reference data. Default column name is
  'reference_names'.

  In table_names, 'reference_version' is a string containing the name of
  the column in 'table' that contains the reference versions. Default
  column name is 'reference_versions'.

  In table_names, 'reference_usage' is a string containing the name of
  the column in 'table' that contains the reference usages. Default
  column name is 'reference_usages'.

  In table_names, 'reference_note' is a string containing the name of
  the column in 'table' that contains the reference notes. Default
  column name is 'reference_notes'.

  In table_names, 'reference_url' is a string containing the name of the
  column in 'table' that contains the reference urls. Default column
  name is 'reference_urls'.

- reference, :

  a list created by the function \[new_reference\]. Optional.

- verbose, :

  boolean indicating whether or not you want progress messages. Default
  = TRUE.

## Value

double - The number of items added

## Examples

``` r
# Create a new empty dataset named 'example_dataset'
data <- new_dataset(dataset_name = "example_dataset")

# Read FASTA data into data.frame
fasta_data <- read_fasta(fasta = strollur_example("final.fasta.gz"))

# Add FASTA sequence data
add(data = data, table = fasta_data, type = "sequences")
#> ℹ Added 2425 sequences.
#> [1] 2425

# To add FASTA data with a resource reference

# Create a new empty dataset named 'example_dataset'
data <- new_dataset(dataset_name = "example_dataset")

# Create a resource reference for the FASTA data
resource_url <- "https://mothur.org/wiki/silva_reference_files/"
resource_reference <-
  new_reference(
    reference_name = "silva.bacteria.fasta",
    reference_version = "1.38.1",
    reference_usage = "alignment by mothur2 v1.0",
    reference_note = "default options",
    reference_url = resource_url
  )

# Add FASTA data with a resource reference

add(
  data = data,
  table = fasta_data,
  type = "sequences",
  reference = resource_reference
)
#> ℹ Added 2425 sequences.
#> [1] 2425

# Add contigs assembly report with a 'sequence_name' column named 'Name'

contigs_report <- readr::read_tsv(
  strollur_example(
    "final.contigs_report.gz"
  ),
  col_names = TRUE, show_col_types = FALSE
)

add(
  data = data, table = contigs_report, type = "reports",
  report_type = "contigs_report", list(sequence_name = "Name")
)
#> ℹ Added a contigs_report.
#> [1] 1

# To add metadata related to your study

metadata <- readr::read_tsv(strollur_example("mouse.dpw.metadata"),
  col_names = TRUE, show_col_types = FALSE
)

add(data = data, table = metadata, type = "metadata")
#> ℹ Added metadata.
#> [1] 1
```
