# Assign sequence abundances, sequence classifications, bins, bin representative sequences, bin classifications, sample distances or treatments to a [strollur object](https://mothur.org/strollur/reference/strollur.html)

Assign sequence abundances, sequence classifications, bins, bin
representative sequences, bin classifications, sample distances or
treatments to a
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Usage

``` r
assign(
  data,
  table,
  type = "bin",
  bin_type = "otu",
  table_names = list(sequence_name = "sequence_name", quality_score = "quality_score",
    abundance = "abundance", sample = "sample", treatment = "treatment", taxonomy =
    "taxonomy", level = "level", confidence = "confidence", bin_name = "bin_name"),
  reference = NULL,
  verbose = TRUE
)
```

## Arguments

- data:

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object

- table:

  a data.frame containing the data you wish to assign

- type:

  a string containing the type of data. Options include:
  `sequence_abundance`, `sequence_taxonomy`, `bin`, `quality`,
  `bin_representative`, `bin_taxonomy`, `sample_distance` and
  `treatment`. Default = `bin`.

- bin_type:

  string containing the bin type you would like the number of bins for.
  Default = `otu`.

- table_names:

  named list used to indicate the names of the columns in the table. By
  default:

  table_names \<- list(sequence_name = "sequence_name", quality_score =
  "quality_score", abundance = "abundance", sample = "sample", treatment
  = "treatment", taxonomy = "taxonomy", level = "level", confidence =
  "confidence", bin_name = "bin_name")

  In table_names, `sequence_name` is a string containing the name of the
  column in 'table' that contains the sequence names. Default column
  name is 'sequence_name'.

  In table_names, `quality_score` is a string containing the name of the
  column in 'table' that contains the sequence quality scores. It is
  used when you are adding quality data. Default column name is
  'quality_score'.

  In table_names, `abundance` is a string containing the name of the
  column in 'table' that contains the abundances. Default column name is
  'abundance'.

  In table_names, `sample` is a string containing the name of the column
  in 'table' that contains the samples. Default column name is 'sample'.

  In table_names, `treatment` is a string containing the name of the
  column in 'table' that contains the treatment names. Default column
  name is 'treatment'.

  In table_names, `taxonomy` is a string containing the name of the
  column in 'table' that contains the classifications. Default column
  name is 'taxonomy'.

  In table_names, `level` is a string containing the name of the column
  in the data.frame that contains the classifications. Default column
  name is `level`.

  In table_names, `confidence` is a string containing the name of the
  column in the data.frame that contains the classifications confidence
  scores. Default column name is `confidence`.

  In table_names, `bin_name` is a string containing the name of the
  column in 'table' that contains the bin names. Default column name is
  'bin_name'.

- reference:

  a list created by the function
  [new_reference](https://mothur.org/strollur/reference/new_reference.md).
  Optional.

- verbose:

  boolean indicating whether or not you want progress messages. Default
  = TRUE.

## Value

an updated
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## See also

[taxonomy_table_format](https://mothur.org/strollur/reference/taxonomy_table_format.md)

## Examples

``` r

# Assign sequence classifications

# create a new empty strollur object named 'example_dataset'
data <- strollur::new_dataset(dataset_name = "example_dataset")

sequence_classifications <- strollur::read_mothur_taxonomy(strollur_example(
  "final.taxonomy.gz"
))

strollur::assign(
  data,
  table = sequence_classifications, type = "sequence_taxonomy"
)
#> Assigned 2425 sequence taxonomies.

# Assigning bins

# read mothur's otu list file into data.frame
otu_data <- strollur::read_mothur_list(list = strollur_example(
  "final.opti_mcc.list.gz"
))

# read mothur's asv list file into data.frame
asv_data <- strollur::read_mothur_list(list = strollur_example(
  "final.asv.list.gz"
))

# read mothur's phylotype list file into data.frame
phylo_data <- strollur::read_mothur_list(list = strollur_example(
  "final.tx.list.gz"
))

# read otu bin representative sequences into a data.frame
bin_reps <- readRDS(strollur_example("miseq_representative_sequences.rds"))

# assign 'otu' bins using sequence names
strollur::assign(data, table = otu_data, bin_type = "otu")
#> Assigned 531 otu bins.

# assign 'asv' bins using sequence names
strollur::assign(data, table = asv_data, bin_type = "asv")
#> Assigned 2425 asv bins.

# assign 'phylotype' bins using sequence names
strollur::assign(data, table = phylo_data, bin_type = "phylotype")
#> Assigned 63 phylotype bins.

# assign 'otu' bin representative sequences
strollur::assign(data, table = bin_reps, type = "bin_representative")
#> Assigned 531 otu bin representative sequences.

# To assign abundance only bins

# create a new empty strollur object named 'example_dataset'
data <- strollur::new_dataset(dataset_name = "example_dataset")

# read mothur's shared file
otu_data <-
  strollur::read_mothur_shared(strollur_example("final.opti_mcc.shared"))

# assign abundance only otus parsed by sample
strollur::assign(data, table = otu_data, bin_type = "otu")
#> Assigned 531 otu bins.

# Assigning bin classifications

# read bin taxonomies
otu_data <- strollur::read_mothur_cons_taxonomy(strollur_example(
  "final.cons.taxonomy"
))

# assign otu consensus taxonomies
strollur::assign(
  data,
  table = otu_data,
  type = "bin_taxonomy", bin_type = "otu"
)
#> Assigned 531 otu bin taxonomies.

# Assign treatments

sample_assignments <- readRDS(strollur_example("miseq_sample_design.rds"))

strollur::assign(data, table = sample_assignments, type = "treatment")
#> Assigned 19 samples to treatments.

# Assign sample distances

dist_file <- strollur_example("final.opti_mcc.jclass.0.03.column.dist")
sample_dists <- readr::read_table(dist_file,
  col_names = FALSE,
  show_col_types = FALSE
)
strollur::assign(data, table = sample_dists, type = "sample_distance")
#> Assigned 171 samples distances.
```
