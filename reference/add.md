# Add sequences, reports, metadata or resource references to a [strollur](https://mothur.org/strollur/reference/strollur.md) object

Add sequences, reports, metadata or resource references to a
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Usage

``` r
add(
  data,
  table,
  type = "sequence",
  report_type = NULL,
  table_names = list(sequence_name = "sequence_name", sequence = "sequence", comment =
    "comment", reference_vendor = "vendor", reference_name = "name", reference_version =
    "version", reference_usage = "usage", reference_note = "note", reference_method_url =
    "method_url", reference_documentation_url = "documentation_url", reference_parameter
    = "parameter", reference_citation = "citation"),
  reference = NULL,
  verbose = TRUE
)
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object

- table, :

  a data.frame containing the data you wish to add.

- type, :

  a string containing the type of data. Options include: 'sequence',
  'resource_reference' 'metadata' and 'report'.

- report_type, :

  a string containing the type of report you are adding. Options
  include: 'metadata' and custom reports.

- table_names, :

  named list used to indicate the names of the columns in the table. By
  default:

  table_names \<- list(sequence_name = "sequence_name", comment =
  "comment", sequence = "sequence", reference_name = "name",
  reference_vendor = "vendor", reference_version = "version",
  reference_usage = "usage", reference_note = "note",
  reference_documentation_url = "documentation_url",
  reference_method_url = "method_url", reference_parameter =
  "parameter", reference_citation = "citation")

  In table_names, 'sequence_name' is a string containing the name of the
  column in 'table' that contains the sequence names. It is used when
  you are adding FASTA data. Default column name is 'sequence_name'.

  In table_names, 'sequence' is a string containing the name of the
  column in 'table' that contains the sequence nucleotide strings. It is
  used when you are adding FASTA data. Default column name is
  'sequence'.

  In table_names, 'comment' is a string containing the name of the
  column in 'table' that contains the sequence comments. It is used when
  you are adding FASTA data. Default column name is 'comment'.

  In table_names, 'reference_vendor' is a string containing the name of
  the column in 'table' that contains the reference vendor names. It is
  used when you are adding reference data. Default column name is
  'vendor'. In table_names, 'reference_name' is a string containing the
  name of the column in 'table' that contains the reference names. It is
  used when you are adding reference data. Default column name is
  'name'.

  In table_names, 'reference_version' is a string containing the name of
  the column in 'table' that contains the reference versions. Default
  column name is 'version'.

  In table_names, 'reference_usage' is a string containing the name of
  the column in 'table' that contains the reference usages. Default
  column name is 'usage'.

  In table_names, 'reference_note' is a string containing the name of
  the column in 'table' that contains the reference notes. Default
  column name is 'note'.

  In table_names, 'reference_method_url' is a string containing the name
  of the column in 'table' that contains the reference method urls.
  Default column name is 'method_url'.

  In table_names, 'reference_documentation_url' is a string containing
  the name of the column in 'table' that contains the reference urls.
  Default column name is 'documentation_url'.

  In table_names, 'reference_parameter' is a string containing the name
  of the column in 'table' that contains the reference parameters.
  Default column name is 'parameter'.

  In table_names, 'reference_citation' is a string containing the name
  of the column in 'table' that contains the reference citations.
  Default column name is 'citation'.

- reference, :

  a list created by the function \[new_reference\]. Optional.

- verbose, :

  boolean indicating whether or not you want progress messages. Default
  = TRUE.

## Value

double - The number of items added

## Examples

``` r

# Create a new empty strollur object named 'example_dataset'
data <- new_dataset(dataset_name = "example_dataset")

# Read FASTA data into data.frame
fasta_data <- read_fasta(fasta = strollur_example("final.fasta.gz"))

# Add FASTA sequence data
add(data = data, table = fasta_data, type = "sequence")
#> Added 2425 sequences.
#> [1] 2425

# To add FASTA data with a resource reference

# Create a new empty strollur object named 'example_dataset'
data <- new_dataset(dataset_name = "example_dataset")

# Create a resource reference for the FASTA data silva_resource <-
silva_resource <- new_reference( vendor = "SILVA", name =
"silva.bacteria.fasta", version = "1.38.1", usage = "alignment of sequences",
note = "reference trimmed to V4 region", method_url =
"https://mothur.org/blog/2024/SILVA-v138_2-reference-files/",
documentation_url = "https://mothur.org/wiki/silva_reference_files/" )

# Add FASTA data with a resource reference

add(
  data,
  table = fasta_data,
  type = "sequence",
  reference = silva_resource
)
#> Added 2425 sequences.
#> Added 1 resource references.
#> [1] 2425

# Add contigs assembly report with a 'sequence_name' column named 'Name'

contigs_report <- readRDS(strollur_example("miseq_contigs_report.rds"))

add(
  data,
  table = contigs_report, type = "report",
  report_type = "contigs_report", list(sequence_name = "Name")
)
#> Added a contigs_report.
#> [1] 1

# To add metadata related to your study

metadata <- readRDS(strollur_example("miseq_metadata.rds"))

add(data, table = metadata, type = "metadata")
#> Added metadata.
#> [1] 1
```
