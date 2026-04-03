# Assign sequence abundances, sequence classifications, bins, bin representative sequences, bin classifications or treatments to a [strollur](https://mothur.org/strollur/reference/strollur.html) object

Assign sequence abundances, sequence classifications, bins, bin
representative sequences, bin classifications or treatments to a
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Usage

``` r
assign(
  data,
  table,
  type = "bins",
  bin_type = "otu",
  table_names = list(sequence_name = "sequence_names", abundance = "abundances", sample =
    "samples", treatment = "treatments", taxonomy = "taxonomies", bin_name = "bin_names"),
  reference = NULL,
  verbose = TRUE
)
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object

- table, :

  a data.frame containing the data you wish to assign

- type, :

  a string containing the type of data. Options include:
  'sequence_abundance', 'sequence_taxonomy', 'bins',
  'bin_representatives', 'bin_taxonomy' and 'treatments'. Default =
  "bins".

- bin_type, :

  string containing the bin type you would like the number of bins for.
  Default = "otu".

- table_names, :

  named list used to indicate the names of the columns in the table. By
  default:

  table_names \<- list(sequence_name = "sequence_names", abundance =
  "abundances", sample = "samples", treatment = "treatments", taxonomy =
  "taxonomies", bin_name = "bin_names")

  In table_names, 'sequence_name' is a string containing the name of the
  column in 'table' that contains the sequence names. Default column
  name is 'sequence_names'.

  In table_names, 'abundance' is a string containing the name of the
  column in 'table' that contains the abundances. Default column name is
  'abundances'.

  In table_names, 'sample' is a string containing the name of the column
  in 'table' that contains the samples. Default column name is
  'samples'.

  In table_names, 'treatment' is a string containing the name of the
  column in 'table' that contains the treatment names. Default column
  name is 'treatments'.

  In table_names, 'taxonomy' is a string containing the name of the
  column in 'table' that contains the classifications. Default column
  name is 'taxonomies'.

  In table_names, 'bin_name' is a string containing the name of the
  column in 'table' that contains the bin names. Default column name is
  'bin_names'.

- reference, :

  a list created by the function \[new_reference\]. Optional.

- verbose, :

  boolean indicating whether or not you want progress messages. Default
  = TRUE.

## Value

double - The number of items assigned

## Examples

``` r
# Assign sequence classifications

# create a new empty strollur object named 'example_dataset'
data <- new_dataset(dataset_name = "example_dataset")

sequence_classifications <- read_mothur_taxonomy(strollur_example(
  "final.taxonomy.gz"
))

assign(
  data,
  table = sequence_classifications, type = "sequence_taxonomy"
)
#> ℹ Assigned 2425 sequence taxonomies.
#> [1] 2425

# Assigning bins

# read mothur's otu list file into data.frame
otu_data <- read_mothur_list(list = strollur_example(
  "final.opti_mcc.list.gz"
))

# read mothur's asv list file into data.frame
asv_data <- read_mothur_list(list = strollur_example(
  "final.asv.list.gz"
))

# read mothur's phylotype list file into data.frame
phylo_data <- read_mothur_list(list = strollur_example(
  "final.tx.list.gz"
))

# read otu bin representative sequences into a data.frame
bin_reps <- readRDS(strollur_example("miseq_representative_sequences.rds"))

# assign 'otu' bins using sequence names
assign(data, table = otu_data, bin_type = "otu")
#> ℹ Assigned 531 otu bins.
#> [1] 531

# assign 'asv' bins using sequence names
assign(data, table = asv_data, bin_type = "asv")
#> ℹ Assigned 2425 asv bins.
#> [1] 2425

# assign 'phylotype' bins using sequence names
assign(data, table = phylo_data, bin_type = "phylotype")
#> ℹ Assigned 63 phylotype bins.
#> [1] 63

# assign 'otu' bin representative sequences
assign(data, table = bin_reps, type = "bin_representatives")
#> ℹ Assigned 531 otu bin representative sequences.
#> [1] 531

# To assign abundance only bins

# create a new empty strollur object named 'example_dataset'
data <- new_dataset(dataset_name = "example_dataset")

# read mothur's shared file
otu_data <- read_mothur_shared(strollur_example("final.opti_mcc.shared"))

# assign abundance only otus parsed by sample
assign(data, table = otu_data, bin_type = "otu")
#> ℹ Assigned 531 otu bins.
#> [1] 531

# Assigning bin classifications

# read bin taxonomies
otu_data <- read_mothur_cons_taxonomy(strollur_example(
  "final.cons.taxonomy"
))

# assign otu consensus taxonomies
assign(
  data,
  table = otu_data,
  type = "bin_taxonomy", bin_type = "otu"
)
#> ℹ Assigned 531 otu bin taxonomies.
#> [1] 531

# Assign treatments

sample_assignments <- readRDS(strollur_example("miseq_sample_design.rds"))

assign(data, table = sample_assignments, type = "treatments")
#> ℹ Assigned 19 samples to treatments.
#> [1] 19
```
