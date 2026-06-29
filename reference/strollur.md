# The \`strollur\` object stores the data associated with your amplicon sequence analysis.

'strollur' is an R6 class that stores nucleotide sequences, abundance,
sample and treatment assignments, taxonomic classifications, asv / otu
clusters and various reports. It is designed to facilitate data analysis
across multiple R packages.

## Author

Sarah Westcott, <swestcot@umich.edu>

## Public fields

- `data`:

  Rcpp::XPtr\<Dataset\> pointer to 'Dataset' c++ class. This allows
  package developers an easy access point to the underlying C++ code
  with additional functionality.

- `raw`:

  Rcpp::RawVector containing the serialized data of the 'Dataset' c++
  class. This allows the load and save functions to work with the class.

- `sequence_tree`:

  a tree that relates sequences to each other

- `sample_tree`:

  a tree that relates samples to each other

## Methods

### Public methods

- [`strollur$new()`](#method-strollur-initialize)

- [`strollur$print()`](#method-strollur-print)

- [`strollur$abundance()`](#method-strollur-abundance)

- [`strollur$add()`](#method-strollur-add)

- [`strollur$add_sample_tree()`](#method-strollur-add_sample_tree)

- [`strollur$add_sequence_tree()`](#method-strollur-add_sequence_tree)

- [`strollur$assign()`](#method-strollur-assign)

- [`strollur$clear()`](#method-strollur-clear)

- [`strollur$count()`](#method-strollur-count)

- [`strollur$get_bin_types()`](#method-strollur-get_bin_types)

- [`strollur$get_sample_tree()`](#method-strollur-get_sample_tree)

- [`strollur$get_sequence_tree()`](#method-strollur-get_sequence_tree)

- [`strollur$get_version()`](#method-strollur-get_version)

- [`strollur$is_equal()`](#method-strollur-is_equal)

- [`strollur$names()`](#method-strollur-names)

- [`strollur$report()`](#method-strollur-report)

- [`strollur$summary()`](#method-strollur-summary)

- [`strollur$clone()`](#method-strollur-clone)

------------------------------------------------------------------------

### `strollur$new()`

Create a new strollur dataset

#### Usage

    strollur$new(name = "", dataset = NULL)

#### Arguments

- `name`:

  String, name of dataset (optional)

- `dataset`:

  a \`strollur\` object.

#### Returns

A new \`strollur\` object.

#### Examples

    # to create an empty strollur object, run the following:

    data <- new_dataset("soil")

------------------------------------------------------------------------

### `strollur$print()`

Print summary of \`strollur\` object

#### Usage

    strollur$print()

#### Returns

No return value, called for side effects.

#### Examples

    miseq <- load_dataset(strollur_example("miseq_sop.rds"))
    miseq

------------------------------------------------------------------------

### `strollur$abundance()`

Get the abundance data for sequences, bins, samples, and treatments.

#### Usage

    strollur$abundance(type = "sequence", bin_type = "otu", by_sample = FALSE)

#### Arguments

- `type`:

  string containing the type of data you want the number of. Options
  include: "sequence", "bin", "sample" and "treatment". Default =
  "sequence".

- `bin_type`:

  string containing the bin type you would like the abundance data for.
  Default = "otu".

- `by_sample`:

  Boolean. When by_sample is TRUE, the abundance data will be parsed by
  sample. Default = FALSE.

#### Returns

data.frame

#### Examples

    miseq <- load_dataset(strollur_example("miseq_sop.rds"))

    # To the total abundance for each sequence
    miseq$abundance(type = "sequence") |> head(n = 5)

    # To the total abundance for each sequence parsed by sample
    miseq$abundance(type = "sequence", by_sample = TRUE) |> head(n = 5)

    # To the total abundance for each "otu" bin
    miseq$abundance(type = "bin", bin_type = "otu") |> head(n = 5)

    # To the total abundance for each "otu" bin parsed by sample
    miseq$abundance(type = "bin", bin_type = "otu", by_sample = TRUE) |>
    head(n = 5)

    # To the total abundance for each "asv" bin
    miseq$abundance(type = "bin", bin_type = "asv") |> head(n = 5)

    # To the total abundance for each "asv" bin parsed by sample
    miseq$abundance(type = "bin", bin_type = "asv", by_sample = TRUE) |>
    head(n = 5)

    # To the total abundance for each sample
    miseq$abundance(type = "sample") |> head(n = 5)

    # To the total abundance for each treatment
    miseq$abundance(type = "treatment")

------------------------------------------------------------------------

### `strollur$add()`

Add sequences, reports or resource references

#### Usage

    strollur$add(
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

#### Arguments

- `table`:

  a data.frame containing the data you wish to add.

- `type`:

  a string containing the type of data. Options include: 'sequence',
  'resource_reference' and 'report'.

- `report_type`:

  a string containing the type of report you are adding.

- `table_names`:

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
  used when ' you are adding reference data. Default column name is
  'vendor'.

  In table_names, 'reference_name' is a string containing the name of
  the ' column in 'table' that contains the reference names. It is used
  when you are ' adding reference data. Default column name is 'name'.

  In table_names, 'reference_version' is a string containing the name of
  the ' column in 'table' that contains the reference versions. Default
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

- `reference`:

  a list created by the function \[new_reference\]. Optional.

- `verbose`:

  boolean indicating whether or not you want progress messages. Default
  = TRUE.

#### Returns

Updated \`strollur\` object - invisible(self)

#### Examples

    fasta_data <- read_fasta(fasta = strollur_example("final.fasta.gz"))
    contigs_report <- readRDS(strollur_example("miseq_contigs_report.rds"))

    # Create a new empty `strollur` object named 'example_dataset'
    data <- new_dataset(dataset_name = "example_dataset")

    data$add(table = fasta_data, type = "sequence")
    data$add(
      table = contigs_report, type = "report",
      report_type = "contigs_report", list(sequence_name = "Name")
    )

    # To add metadata related to your study

    metadata <- readRDS(strollur_example("miseq_metadata.rds"))

    data$add(table = metadata, type = "report", report_type = "metadata")

------------------------------------------------------------------------

### `strollur$add_sample_tree()`

Add phylo tree relating the samples in your dataset

#### Usage

    strollur$add_sample_tree(tree)

#### Arguments

- `tree`:

  a phylo tree object created by ape::read.tree.

#### Returns

Updated \`strollur\` object

#### Examples

     data <- new_dataset("my_dataset")

     df <- read_mothur_shared(strollur_example("final.opti_mcc.shared"))
     assign(data = data, table = df, type = "bin", bin_type = "otu")

     tree <- ape::read.tree(strollur_example(
     "final.opti_mcc.jclass.ave.tre"))

     data$add_sample_tree(tree)

------------------------------------------------------------------------

### `strollur$add_sequence_tree()`

Add phylo tree relating the sequences in your dataset

#### Usage

    strollur$add_sequence_tree(tree)

#### Arguments

- `tree`:

  a phylo tree object created by ape::read.tree.

#### Returns

Updated \`strollur\` object

#### Examples

     data <- new_dataset("my_dataset")
     tree <- ape::read.tree(strollur_example("final.phylip.tre.gz"))
     data$add_sequence_tree(tree)

------------------------------------------------------------------------

### `strollur$assign()`

Assign sequence abundances, sequence classifications, bins, bin
representative sequences, bin classifications or treatments.

#### Usage

    strollur$assign(
      table,
      type = "bin",
      bin_type = "otu",
      table_names = list(sequence_name = "sequence_name", abundance = "abundance", sample =
        "sample", treatment = "treatment", taxonomy = "taxonomy", bin_name = "bin_name"),
      reference = NULL,
      verbose = TRUE
    )

#### Arguments

- `table`:

  a data.frame containing the data you wish to assign

- `type`:

  a string containing the type of data. Options include:
  'sequence_abundance', 'sequence_taxonomy', 'bin',
  'bin_representative', 'bin_taxonomy' and 'treatment'. Default = "bin".

- `bin_type`:

  string containing the bin type you would like the number of bins for.
  Default = "otu".

- `table_names`:

  named list used to indicate the names of the columns in the table. By
  default:

  table_names \<- list(sequence_name = "sequence_name", abundance =
  "abundance", sample = "sample", treatment = "treatment", taxonomy =
  "taxonomy", bin_name = "bin_name")

  In table_names, 'sequence_name' is a string containing the name of the
  column in 'table' that contains the sequence names. Default column
  name is 'sequence_name'.

  In table_names, 'abundance' is a string containing the name of the
  column in 'table' that contains the abundances. Default column name is
  'abundance'.

  In table_names, 'sample' is a string containing the name of the column
  in 'table' that contains the samples. Default column name is 'sample'.

  In table_names, 'treatment' is a string containing the name of the
  column in 'table' that contains the treatment names. Default column
  name is 'treatment'.

  In table_names, 'taxonomy' is a string containing the name of the
  column in 'table' that contains the classifications. Default column
  name is 'taxonomy'.

  In table_names, 'bin_name' is a string containing the name of the
  column in 'table' that contains the bin names. Default column name is
  'bin_name'.

- `reference`:

  a list created by the function \[new_reference\]. Optional.

- `verbose`:

  boolean indicating whether or not you want progress messages. Default
  = TRUE.

#### Returns

Updated \`strollur\` object

#### Examples

    # create a new empty strollur object named 'example_dataset'

    data <- new_dataset(dataset_name = "example_dataset")

    # Assign sequence abundances

    abundance_by_sample <- read_mothur_count(strollur_example(
      "final.count_table.gz"
    ))

    data$assign(table = abundance_by_sample, type = "sequence_abundance")

    # Assign sequence classifications

    sequence_classifications <- read_mothur_taxonomy(strollur_example(
      "final.taxonomy.gz"
    ))

    data$assign(table = sequence_classifications, type = "sequence_taxonomy")

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
    bin_reps <- readRDS(strollur_example(
                            "miseq_representative_sequences.rds"))

    # assign 'otu' bins using sequence names
    data$assign(table = otu_data, bin_type = "otu")

    # assign 'asv' bins using sequence names
    data$assign(table = asv_data, bin_type = "asv")

    # assign 'phylotype' bins using sequence names
    data$assign(table = phylo_data, bin_type = "phylotype")

    # assign 'otu' bin representative sequences
    data$assign(table = bin_reps, type = "bin_representative")

    # To assign abundance only bins

    # create a new empty strollur object named 'example_dataset'
    data <- new_dataset(dataset_name = "example_dataset")

    # read mothur's shared file
    otu_data <- read_mothur_shared(strollur_example("final.opti_mcc.shared"))

    # assign abundance only otus parsed by sample
    data$assign(table = otu_data, bin_type = "otu")

    # Assigning bin classifications

    # read bin taxonomies
    otu_data <- read_mothur_cons_taxonomy(strollur_example(
      "final.cons.taxonomy"
    ))

    # assign otu consensus taxonomies
    data$assign(
      table = otu_data,
      type = "bin_taxonomy", bin_type = "otu"
    )

    # Assign treatments

    sample_assignments <- readRDS(
       strollur_example("miseq_sample_design.rds"))

    data$assign(table = sample_assignments, type = "treatment")

------------------------------------------------------------------------

### `strollur$clear()`

Clear data from datasest

#### Usage

    strollur$clear()

#### Returns

Updated \`strollur\` object

#### Examples

    miseq <- load_dataset(strollur_example("miseq_sop.rds"))
    miseq
    miseq$clear()
    miseq

------------------------------------------------------------------------

### `strollur$count()`

Find the number of sequences, samples, treatments or bins of a given
type

#### Usage

    strollur$count(
      type = "sequence",
      bin_type = "otu",
      samples = NULL,
      distinct = FALSE
    )

#### Arguments

- `type`:

  string containing the type of data you want the number of. Options
  include: "sequence", "sample", "treatment", "bin", and
  "resource_reference". Default = "sequence".

- `bin_type`:

  string containing the bin type you would like the number of bins for.
  Default = "otu".

- `samples`:

  vector of strings. samples is only used when 'type' = "sequence" or
  'type' = "bin" . samples should contain the names of the samples you
  want the count for. Default = NULL.

- `distinct`:

  Boolean. distinct is used when 'type' = "sequence" or 'type' = "bin".
  When 'type' = "sequence" and distinct is TRUE the number of unique
  sequences is returned. When 'type' = "sequence" and distinct is FALSE
  the total number of sequences is returned. This can also be combined
  with samples to find the number of unique sequences found ONLY in a
  given set of samples, or to find the number of unique sequences in
  given set of samples that may also be present in other samples. When
  'type' = "bin", you can set distinct = TRUE to return the number of
  bins that ONLY contain sequences from the given samples. When distinct
  is FALSE the count returned contains bins with sequences from a given
  samples, but those bins may also contain other samples. Default =
  FALSE.

#### Returns

double

#### Examples

    miseq <- load_dataset(strollur_example("miseq_sop.rds"))

    # To get the total number of sequences
    miseq$count(type = "sequence")

    # To get number of unique sequences
    miseq$count(type = "sequence", distinct = TRUE)

    # To get number of unique sequences from samples 'F3D0' and 'F3D1'
    # Note these sequences will be present in both samples but may be
    # be present in other samples as well
    miseq$count(type = "sequence", samples = c("F3D0", "F3D1"))

    # To get number of unique sequences exclusive to samples 'F3D0' and
    # 'F3D1'. Note sequences are present in both samples and NOT present in
    # any other samples.

    miseq$count(type = "sequence",
                samples = c("F3D0", "F3D1"), distinct = TRUE )

    # To get the number of samples in the dataset
    miseq$count(type = "sample")

    # To get the number of treatments in the dataset
    miseq$count(type = "treatment")

    # To get the number of "otu" bins in the dataset
    miseq$count(type = "bin", bin_type = "otu")

    # To get the number of "asv" bins in the dataset
    miseq$count(type = "bin", bin_type = "asv")

    # To get the number of "phylotype" bins in the dataset
    miseq$count(type = "bin", bin_type = "phylotype")

    # To get number of "otu" bins from samples 'F3D0' and 'F3D1'
    # Note these bins will have sequences from both samples but there may be
    # other samples present as well
    miseq$count(
      type = "bin", bin_type = "otu", samples = c("F3D0", "F3D1")
    )

    # To get number of "otu" bins unique to samples 'F3D0' and 'F3D1'
    # Note these bins will have sequences from both samples and NO other
    # samples will be present in the bins.

    miseq$count(
      type = "bin", bin_type = "otu",
      samples = c("F3D0", "F3D1"), distinct = TRUE
    )

------------------------------------------------------------------------

### `strollur$get_bin_types()`

Get bin table types

#### Usage

    strollur$get_bin_types()

#### Returns

vector of strings

#### Examples

    data <- miseq_sop_example()
    data$get_bin_types()

------------------------------------------------------------------------

### `strollur$get_sample_tree()`

Get phylo tree relating the samples in your dataset.

#### Usage

    strollur$get_sample_tree()

#### Returns

ape::tree

#### Examples

     tree <- ape::read.tree(strollur_example(
      "final.opti_mcc.jclass.ave.tre"))

     df <- read_mothur_shared(strollur_example("final.opti_mcc.shared"))

     data <- new_dataset("my_dataset")

     # assign abundance 'otu' bins
     data$assign(table = df, type = "bin", bin_type = "otu")

     data$add_sample_tree(tree)
     data$get_sample_tree()

------------------------------------------------------------------------

### `strollur$get_sequence_tree()`

Get phylo tree relating the sequences in your strollur object.

#### Usage

    strollur$get_sequence_tree()

#### Returns

ape::tree

#### Examples

     data <- new_dataset("my_dataset")
     tree <- ape::read.tree(strollur_example("final.phylip.tre.gz"))
     data$add_sequence_tree(tree)
     data$get_sequence_tree()

------------------------------------------------------------------------

### `strollur$get_version()`

Get the version of the
[strollur](https://mothur.org/strollur/reference/strollur.html) object.

#### Usage

    strollur$get_version()

#### Returns

a logical

#### Examples

    data <- new_dataset("test")

    data$get_version()

------------------------------------------------------------------------

### `strollur$is_equal()`

Determine if two
[strollur](https://mothur.org/strollur/reference/strollur.html) objects
are equal.

#### Usage

    strollur$is_equal(data)

#### Arguments

- `data, `:

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object

#### Returns

a logical

#### Examples

    miseq <- load_dataset(strollur_example("miseq_sop.rds"))

    data <- copy_dataset(miseq)

    miseq$is_equal(data)

------------------------------------------------------------------------

### `strollur$names()`

Get the names of a given type of data

#### Usage

    strollur$names(
      type = "sequence",
      bin_type = "otu",
      samples = NULL,
      distinct = FALSE
    )

#### Arguments

- `type`:

  string containing the type of data you would like. Options include:
  "dataset", "sequence", "bin", "sample", "treatment", "report". Default
  = "sequence".

- `bin_type`:

  string containing the bin type you would like the names for. Default =
  "otu".

- `samples`:

  vector of strings. samples is only used when 'type' = "sequence" or
  'type' = "bin" . samples should contain the names of the samples you
  want names for. Default = NULL.

- `distinct`:

  Boolean. distinct is used when 'type' = "sequence" or 'type' = "bin"
  and the samples parameter is used. The distinct parameter allows you
  to get the names that present given set of samples. When distinct is
  TRUE, the names function will return the names that ONLY contain data
  from the given samples. When distinct is FALSE the data returned
  contains data from a given samples, but may ALSO contain data from
  other samples. Default = FALSE.

#### Returns

vector of strings, containing the names requested

#### Examples

    miseq <- load_dataset(strollur_example("miseq_sop.rds"))

    # To get the name of the dataset
    miseq$names(type = "dataset")

    # To get the names of the sequences
    miseq$names(type = "sequence")

    # To get the names of the sequences present sample 'F3D0'
    miseq$names(type = "sequence", samples = c("F3D0"))

    #' # To get the names of the sequences unique to sample 'F3D0'
    miseq$names(type = "sequence", samples = c("F3D0"), distinct = TRUE)

    # To get the names of the samples
    miseq$names(type = "sample")

    # To get the names of the treatments
    miseq$names(type = "treatment")

    # To get the names of the bins
    miseq$names(type = "bin")

    # To get the names of the bins that are unique to 'F3D0'
    miseq$names(type = "bin", samples = c("F3D0"), distinct = TRUE)

    # To get the names of the bins that include sequences from 'F3D0'
    miseq$names(type = "bin", samples = c("F3D0"), distinct = FALSE)

    # To get the names of the reports
    miseq$names(type = "report")

------------------------------------------------------------------------

### `strollur$report()`

Get a data.frame containing the given report

#### Usage

    strollur$report(type = "sequence", bin_type = "otu")

#### Arguments

- `type`:

  string containing the type of report you would like. Options include:
  "fasta", "sequence", "sequence_bin_assignment", "sequence_taxonomy",
  "bin_taxonomy", "bin_representative", "sample_assignment",
  "resource_reference", "sequence_scrap", "bin_scrap". If you have added
  custom reports for alignment, contigs_assembly or chimeras, you can
  get those as well. Default = "sequence".

- `bin_type`:

  string containing the bin type you would like a bin_taxonomy report
  for. Default = "otu".

#### Returns

data.frame

#### Examples

    miseq <- load_dataset(strollur_example("miseq_sop.rds"))

    # To get the FASTA data

    miseq$report(type = "fasta") |> head(n = 5)

    # To get a report about the FASTA data

    miseq$report(type = "sequence") |> head(n = 5)

    # To get the sequence bin assignments

    miseq$report(type = "sequence_bin_assignment", bin_type = "otu") |>
    head(n = 5)

    # To get the sample treatment assignments

    miseq$report(type = "sample_assignment") |> head(n = 5)

    # To get a report about sequence classifications

    miseq$report(type = "sequence_taxonomy") |> head(n = 5)

    # To get a report about bin classifications for 'otu' data

    miseq$report(type = "bin_taxonomy", bin_type = "otu") |> head(n = 5)

    # To get the 'otu' bin representative sequences

    miseq$report(type = "bin_representative", bin_type = "otu") |>
    head(n = 5)

    # To get a report about the sequences removed during your analysis:

    miseq$report(type = "sequence_scrap")

    # To get a report about the "otu" bins removed during your analysis:

    miseq$report(type = "bin_scrap", bin_type = "otu")

    # To get the metadata associated with your data:

    metadata <- miseq$report(type = "metadata") |> head(n = 5)

    # To get the resource references associated with your data:

    references <- miseq$report(type = "resource_reference")

    # To get our custom report containing the contigs assembly data:

    miseq$report(type = "contigs_report") |> head(n = 5)

------------------------------------------------------------------------

### `strollur$summary()`

Summarize the sequences data, custom reports, and scrapped data

#### Usage

    strollur$summary(type = "sequence", report_type = NULL, verbose = TRUE)

#### Arguments

- `type`:

  string containing the type of data you want the number of. Options
  include: "sequence", "report" and "scrap". Default = "sequence".

- `report_type`:

  string containing the report type you would summarized. For example,
  the miseq_sop_example includes contigs assembly data and can be
  accessed with report_type = "contigs_report". Default = NULL.

- `verbose`:

  boolean indicating whether or not you want progress messages. Default
  = TRUE.

#### Returns

data.frame Get summary of the sequence reports

#### Examples

    miseq <- load_dataset(strollur_example("miseq_sop.rds"))

    # To get the summary of your FASTA data
    miseq$summary(type = "sequence")

    # summarize contigs_report
    miseq$summary(type = "report", report_type = "contigs_report")

    # remove sample 'F3D0' to produce a scrap report
    xdev_remove_samples(data = miseq, samples = c("F3D0"))

    # summarize scrapped data -
    # sequences and bins scrapped by removing the sample "F3D0"
    miseq$summary(type = "scrap")

------------------------------------------------------------------------

### `strollur$clone()`

The objects of this class are cloneable with this method.

#### Usage

    strollur$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r

## ------------------------------------------------
## Method `strollur$new()`
## ------------------------------------------------


# to create an empty strollur object, run the following:

data <- new_dataset("soil")


## ------------------------------------------------
## Method `strollur$print()`
## ------------------------------------------------

miseq <- load_dataset(strollur_example("miseq_sop.rds"))
miseq
#> miseq_sop:
#> 
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  375    249      0        3     0      1.00
#> 2.5%-tile:       1  375    252      0        4     0   2849.08
#> 25%-tile:        1  375    252      0        4     0  28490.75
#> Median:          1  375    253      0        4     0  56981.50
#> 75%-tile:        1  375    253      0        5     0  85472.25
#> 97.5%-tile:      1  375    254      0        6     0 111113.93
#> Maximum:         1  375    256      0        6     0 113963.00
#> Mean:            1  375    252      0        4     0  56981.64
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
#> Total number of custom reports: 2 
#> 


## ------------------------------------------------
## Method `strollur$abundance()`
## ------------------------------------------------


miseq <- load_dataset(strollur_example("miseq_sop.rds"))

# To the total abundance for each sequence
miseq$abundance(type = "sequence") |> head(n = 5)
#>                                  sequence_name abundance
#> 1  M00967_43_000000000-A3JHG_1_1101_10133_8460       620
#> 2 M00967_43_000000000-A3JHG_1_1101_10331_23332         1
#> 3 M00967_43_000000000-A3JHG_1_1101_10382_22128         1
#> 4 M00967_43_000000000-A3JHG_1_1101_11035_15765         1
#> 5 M00967_43_000000000-A3JHG_1_1101_11348_22601         2

# To the total abundance for each sequence parsed by sample
miseq$abundance(type = "sequence", by_sample = TRUE) |> head(n = 5)
#>                                 sequence_name abundance sample treatment
#> 1 M00967_43_000000000-A3JHG_1_1101_10133_8460        32   F3D0     Early
#> 2 M00967_43_000000000-A3JHG_1_1101_10133_8460       127   F3D1     Early
#> 3 M00967_43_000000000-A3JHG_1_1101_10133_8460         1 F3D146      Late
#> 4 M00967_43_000000000-A3JHG_1_1101_10133_8460         1 F3D149      Late
#> 5 M00967_43_000000000-A3JHG_1_1101_10133_8460         1 F3D150      Late

# To the total abundance for each "otu" bin
miseq$abundance(type = "bin", bin_type = "otu") |> head(n = 5)
#>   otu_id abundance
#> 1 Otu001     12288
#> 2 Otu002      8892
#> 3 Otu003      7794
#> 4 Otu004      7476
#> 5 Otu005      7450

# To the total abundance for each "otu" bin parsed by sample
miseq$abundance(type = "bin", bin_type = "otu", by_sample = TRUE) |>
head(n = 5)
#>   bin_name abundance sample treatment
#> 1   Otu001       499   F3D0     Early
#> 2   Otu001       351   F3D1     Early
#> 3   Otu001       388 F3D141      Late
#> 4   Otu001       244 F3D142      Late
#> 5   Otu001       189 F3D143      Late

# To the total abundance for each "asv" bin
miseq$abundance(type = "bin", bin_type = "asv") |> head(n = 5)
#>    asv_id abundance
#> 1 Asv0001     12196
#> 2 Asv0002      8829
#> 3 Asv0003      7698
#> 4 Asv0004      7436
#> 5 Asv0005      7307

# To the total abundance for each "asv" bin parsed by sample
miseq$abundance(type = "bin", bin_type = "asv", by_sample = TRUE) |>
head(n = 5)
#>   bin_name abundance sample treatment
#> 1  Asv0001       495   F3D0     Early
#> 2  Asv0001       340   F3D1     Early
#> 3  Asv0001       386 F3D141      Late
#> 4  Asv0001       242 F3D142      Late
#> 5  Asv0001       188 F3D143      Late

# To the total abundance for each sample
miseq$abundance(type = "sample") |> head(n = 5)
#>   sample abundance
#> 1   F3D0      6191
#> 2   F3D1      4652
#> 3 F3D141      4656
#> 4 F3D142      2423
#> 5 F3D143      2403

# To the total abundance for each treatment
miseq$abundance(type = "treatment")
#>   treatment abundance
#> 1     Early     55634
#> 2      Late     58329


## ------------------------------------------------
## Method `strollur$add()`
## ------------------------------------------------


fasta_data <- read_fasta(fasta = strollur_example("final.fasta.gz"))
contigs_report <- readRDS(strollur_example("miseq_contigs_report.rds"))

# Create a new empty `strollur` object named 'example_dataset'
data <- new_dataset(dataset_name = "example_dataset")

data$add(table = fasta_data, type = "sequence")
#> Added 2425 sequences.
data$add(
  table = contigs_report, type = "report",
  report_type = "contigs_report", list(sequence_name = "Name")
)
#> Added a contigs_report report.

# To add metadata related to your study

metadata <- readRDS(strollur_example("miseq_metadata.rds"))

data$add(table = metadata, type = "report", report_type = "metadata")
#> Added a metadata report.


## ------------------------------------------------
## Method `strollur$add_sample_tree()`
## ------------------------------------------------


 data <- new_dataset("my_dataset")

 df <- read_mothur_shared(strollur_example("final.opti_mcc.shared"))
 assign(data = data, table = df, type = "bin", bin_type = "otu")
#> Assigned 531 otu bins.

 tree <- ape::read.tree(strollur_example(
 "final.opti_mcc.jclass.ave.tre"))

 data$add_sample_tree(tree)

## ------------------------------------------------
## Method `strollur$add_sequence_tree()`
## ------------------------------------------------


 data <- new_dataset("my_dataset")
 tree <- ape::read.tree(strollur_example("final.phylip.tre.gz"))
 data$add_sequence_tree(tree)
#> Added 2425 sequences.


## ------------------------------------------------
## Method `strollur$assign()`
## ------------------------------------------------


# create a new empty strollur object named 'example_dataset'

data <- new_dataset(dataset_name = "example_dataset")

# Assign sequence abundances

abundance_by_sample <- read_mothur_count(strollur_example(
  "final.count_table.gz"
))

data$assign(table = abundance_by_sample, type = "sequence_abundance")
#> Assigned 2425 sequence abundances.

# Assign sequence classifications

sequence_classifications <- read_mothur_taxonomy(strollur_example(
  "final.taxonomy.gz"
))

data$assign(table = sequence_classifications, type = "sequence_taxonomy")
#> Assigned 2425 sequence taxonomies.

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
bin_reps <- readRDS(strollur_example(
                        "miseq_representative_sequences.rds"))

# assign 'otu' bins using sequence names
data$assign(table = otu_data, bin_type = "otu")
#> Assigned 531 otu bins.

# assign 'asv' bins using sequence names
data$assign(table = asv_data, bin_type = "asv")
#> Assigned 2425 asv bins.

# assign 'phylotype' bins using sequence names
data$assign(table = phylo_data, bin_type = "phylotype")
#> Assigned 63 phylotype bins.

# assign 'otu' bin representative sequences
data$assign(table = bin_reps, type = "bin_representative")
#> Assigned 531 otu bin representative sequences.

# To assign abundance only bins

# create a new empty strollur object named 'example_dataset'
data <- new_dataset(dataset_name = "example_dataset")

# read mothur's shared file
otu_data <- read_mothur_shared(strollur_example("final.opti_mcc.shared"))

# assign abundance only otus parsed by sample
data$assign(table = otu_data, bin_type = "otu")
#> Assigned 531 otu bins.

# Assigning bin classifications

# read bin taxonomies
otu_data <- read_mothur_cons_taxonomy(strollur_example(
  "final.cons.taxonomy"
))

# assign otu consensus taxonomies
data$assign(
  table = otu_data,
  type = "bin_taxonomy", bin_type = "otu"
)
#> Assigned 531 otu bin taxonomies.

# Assign treatments

sample_assignments <- readRDS(
   strollur_example("miseq_sample_design.rds"))

data$assign(table = sample_assignments, type = "treatment")
#> Assigned 19 samples to treatments.


## ------------------------------------------------
## Method `strollur$clear()`
## ------------------------------------------------

miseq <- load_dataset(strollur_example("miseq_sop.rds"))
miseq
#> miseq_sop:
#> 
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  375    249      0        3     0      1.00
#> 2.5%-tile:       1  375    252      0        4     0   2849.08
#> 25%-tile:        1  375    252      0        4     0  28490.75
#> Median:          1  375    253      0        4     0  56981.50
#> 75%-tile:        1  375    253      0        5     0  85472.25
#> 97.5%-tile:      1  375    254      0        6     0 111113.93
#> Maximum:         1  375    256      0        6     0 113963.00
#> Mean:            1  375    252      0        4     0  56981.64
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
#> Total number of custom reports: 2 
#> 
miseq$clear()
miseq
#> 
#> Total number of seqs: 0 
#> 
#> 

## ------------------------------------------------
## Method `strollur$count()`
## ------------------------------------------------


miseq <- load_dataset(strollur_example("miseq_sop.rds"))

# To get the total number of sequences
miseq$count(type = "sequence")
#> [1] 113963

# To get number of unique sequences
miseq$count(type = "sequence", distinct = TRUE)
#> [1] 2425

# To get number of unique sequences from samples 'F3D0' and 'F3D1'
# Note these sequences will be present in both samples but may be
# be present in other samples as well
miseq$count(type = "sequence", samples = c("F3D0", "F3D1"))
#> [1] 9385

# To get number of unique sequences exclusive to samples 'F3D0' and
# 'F3D1'. Note sequences are present in both samples and NOT present in
# any other samples.

miseq$count(type = "sequence",
            samples = c("F3D0", "F3D1"), distinct = TRUE )
#> [1] 2

# To get the number of samples in the dataset
miseq$count(type = "sample")
#> [1] 19

# To get the number of treatments in the dataset
miseq$count(type = "treatment")
#> [1] 2

# To get the number of "otu" bins in the dataset
miseq$count(type = "bin", bin_type = "otu")
#> [1] 531

# To get the number of "asv" bins in the dataset
miseq$count(type = "bin", bin_type = "asv")
#> [1] 2425

# To get the number of "phylotype" bins in the dataset
miseq$count(type = "bin", bin_type = "phylotype")
#> [1] 63

# To get number of "otu" bins from samples 'F3D0' and 'F3D1'
# Note these bins will have sequences from both samples but there may be
# other samples present as well
miseq$count(
  type = "bin", bin_type = "otu", samples = c("F3D0", "F3D1")
)
#> [1] 125

# To get number of "otu" bins unique to samples 'F3D0' and 'F3D1'
# Note these bins will have sequences from both samples and NO other
# samples will be present in the bins.

miseq$count(
  type = "bin", bin_type = "otu",
  samples = c("F3D0", "F3D1"), distinct = TRUE
)
#> [1] 1


## ------------------------------------------------
## Method `strollur$get_bin_types()`
## ------------------------------------------------


data <- miseq_sop_example()
#> Added 2425 sequences.
#> Assigned 2425 sequence abundances.
#> Assigned 2425 sequence taxonomies.
#> Assigned 531 otu bins.
#> Assigned 2425 asv bins.
#> Assigned 63 phylotype bins.
#> Assigned 19 samples to treatments.
#> Assigned 531 otu bin taxonomies.
#> Assigned 531 otu bin representative sequences.
#> Added a metadata report.
#> Added 2 resource references.
#> Added a contigs_report report.
data$get_bin_types()
#> [1] "otu"       "asv"       "phylotype"


## ------------------------------------------------
## Method `strollur$get_sample_tree()`
## ------------------------------------------------


 tree <- ape::read.tree(strollur_example(
  "final.opti_mcc.jclass.ave.tre"))

 df <- read_mothur_shared(strollur_example("final.opti_mcc.shared"))

 data <- new_dataset("my_dataset")

 # assign abundance 'otu' bins
 data$assign(table = df, type = "bin", bin_type = "otu")
#> Assigned 531 otu bins.

 data$add_sample_tree(tree)
 data$get_sample_tree()
#> 
#> Phylogenetic tree with 19 tips and 18 internal nodes.
#> 
#> Tip labels:
#>   F3D9, F3D8, F3D6, F3D5, F3D2, F3D1, ...
#> 
#> Rooted; includes branch length(s).


## ------------------------------------------------
## Method `strollur$get_sequence_tree()`
## ------------------------------------------------


 data <- new_dataset("my_dataset")
 tree <- ape::read.tree(strollur_example("final.phylip.tre.gz"))
 data$add_sequence_tree(tree)
#> Added 2425 sequences.
 data$get_sequence_tree()
#> 
#> Phylogenetic tree with 2425 tips and 2424 internal nodes.
#> 
#> Tip labels:
#>   M00967_43_000000000-A3JHG_1_1114_15727_25995, M00967_43_000000000-A3JHG_1_2109_19976_22044, M00967_43_000000000-A3JHG_1_1102_9244_9305, M00967_43_000000000-A3JHG_1_2101_14159_9619, M00967_43_000000000-A3JHG_1_1111_12315_7486, M00967_43_000000000-A3JHG_1_1107_12586_26826, ...
#> 
#> Rooted; includes branch length(s).


## ------------------------------------------------
## Method `strollur$get_version()`
## ------------------------------------------------


data <- new_dataset("test")

data$get_version()
#> [1] "0.1.1"


## ------------------------------------------------
## Method `strollur$is_equal()`
## ------------------------------------------------


miseq <- load_dataset(strollur_example("miseq_sop.rds"))

data <- copy_dataset(miseq)

miseq$is_equal(data)
#> → The strollur objects public field 'raw' are not  equivalent.
#> [1] FALSE


## ------------------------------------------------
## Method `strollur$names()`
## ------------------------------------------------


miseq <- load_dataset(strollur_example("miseq_sop.rds"))

# To get the name of the dataset
miseq$names(type = "dataset")
#> [1] "miseq_sop"

# To get the names of the sequences
miseq$names(type = "sequence")
#>    [1] "M00967_43_000000000-A3JHG_1_1101_10133_8460" 
#>    [2] "M00967_43_000000000-A3JHG_1_1101_10331_23332"
#>    [3] "M00967_43_000000000-A3JHG_1_1101_10382_22128"
#>    [4] "M00967_43_000000000-A3JHG_1_1101_11035_15765"
#>    [5] "M00967_43_000000000-A3JHG_1_1101_11348_22601"
#>    [6] "M00967_43_000000000-A3JHG_1_1101_11371_20625"
#>    [7] "M00967_43_000000000-A3JHG_1_1101_11480_8830" 
#>    [8] "M00967_43_000000000-A3JHG_1_1101_11754_11390"
#>    [9] "M00967_43_000000000-A3JHG_1_1101_11983_2593" 
#>   [10] "M00967_43_000000000-A3JHG_1_1101_12302_23776"
#>   [11] "M00967_43_000000000-A3JHG_1_1101_12446_20708"
#>   [12] "M00967_43_000000000-A3JHG_1_1101_12732_20169"
#>   [13] "M00967_43_000000000-A3JHG_1_1101_12861_15386"
#>   [14] "M00967_43_000000000-A3JHG_1_1101_13158_12636"
#>   [15] "M00967_43_000000000-A3JHG_1_1101_13609_10025"
#>   [16] "M00967_43_000000000-A3JHG_1_1101_13714_14338"
#>   [17] "M00967_43_000000000-A3JHG_1_1101_14159_27572"
#>   [18] "M00967_43_000000000-A3JHG_1_1101_14364_8401" 
#>   [19] "M00967_43_000000000-A3JHG_1_1101_14915_26251"
#>   [20] "M00967_43_000000000-A3JHG_1_1101_15011_15369"
#>   [21] "M00967_43_000000000-A3JHG_1_1101_15463_5782" 
#>   [22] "M00967_43_000000000-A3JHG_1_1101_15533_5293" 
#>   [23] "M00967_43_000000000-A3JHG_1_1101_15591_4696" 
#>   [24] "M00967_43_000000000-A3JHG_1_1101_15683_17782"
#>   [25] "M00967_43_000000000-A3JHG_1_1101_16007_27282"
#>   [26] "M00967_43_000000000-A3JHG_1_1101_16143_23069"
#>   [27] "M00967_43_000000000-A3JHG_1_1101_16182_25089"
#>   [28] "M00967_43_000000000-A3JHG_1_1101_16381_21558"
#>   [29] "M00967_43_000000000-A3JHG_1_1101_16390_9529" 
#>   [30] "M00967_43_000000000-A3JHG_1_1101_16914_26242"
#>   [31] "M00967_43_000000000-A3JHG_1_1101_17961_23432"
#>   [32] "M00967_43_000000000-A3JHG_1_1101_18044_1900" 
#>   [33] "M00967_43_000000000-A3JHG_1_1101_18106_21085"
#>   [34] "M00967_43_000000000-A3JHG_1_1101_18143_13375"
#>   [35] "M00967_43_000000000-A3JHG_1_1101_18278_3345" 
#>   [36] "M00967_43_000000000-A3JHG_1_1101_18346_24737"
#>   [37] "M00967_43_000000000-A3JHG_1_1101_18517_10721"
#>   [38] "M00967_43_000000000-A3JHG_1_1101_18682_27059"
#>   [39] "M00967_43_000000000-A3JHG_1_1101_18693_22396"
#>   [40] "M00967_43_000000000-A3JHG_1_1101_18922_4934" 
#>   [41] "M00967_43_000000000-A3JHG_1_1101_19181_15979"
#>   [42] "M00967_43_000000000-A3JHG_1_1101_19534_17052"
#>   [43] "M00967_43_000000000-A3JHG_1_1101_20035_18358"
#>   [44] "M00967_43_000000000-A3JHG_1_1101_2004_15863" 
#>   [45] "M00967_43_000000000-A3JHG_1_1101_20262_22075"
#>   [46] "M00967_43_000000000-A3JHG_1_1101_20589_8650" 
#>   [47] "M00967_43_000000000-A3JHG_1_1101_20702_8361" 
#>   [48] "M00967_43_000000000-A3JHG_1_1101_20815_21052"
#>   [49] "M00967_43_000000000-A3JHG_1_1101_21034_7357" 
#>   [50] "M00967_43_000000000-A3JHG_1_1101_21305_8478" 
#>   [51] "M00967_43_000000000-A3JHG_1_1101_21616_8560" 
#>   [52] "M00967_43_000000000-A3JHG_1_1101_22681_5598" 
#>   [53] "M00967_43_000000000-A3JHG_1_1101_23090_14223"
#>   [54] "M00967_43_000000000-A3JHG_1_1101_23238_24359"
#>   [55] "M00967_43_000000000-A3JHG_1_1101_23714_7638" 
#>   [56] "M00967_43_000000000-A3JHG_1_1101_23774_12021"
#>   [57] "M00967_43_000000000-A3JHG_1_1101_24095_14123"
#>   [58] "M00967_43_000000000-A3JHG_1_1101_24369_12902"
#>   [59] "M00967_43_000000000-A3JHG_1_1101_24693_20989"
#>   [60] "M00967_43_000000000-A3JHG_1_1101_24952_9272" 
#>   [61] "M00967_43_000000000-A3JHG_1_1101_25551_18277"
#>   [62] "M00967_43_000000000-A3JHG_1_1101_25615_24528"
#>   [63] "M00967_43_000000000-A3JHG_1_1101_25878_10062"
#>   [64] "M00967_43_000000000-A3JHG_1_1101_26014_7542" 
#>   [65] "M00967_43_000000000-A3JHG_1_1101_26020_8626" 
#>   [66] "M00967_43_000000000-A3JHG_1_1101_26204_21938"
#>   [67] "M00967_43_000000000-A3JHG_1_1101_26634_17611"
#>   [68] "M00967_43_000000000-A3JHG_1_1101_27744_12678"
#>   [69] "M00967_43_000000000-A3JHG_1_1101_28416_18742"
#>   [70] "M00967_43_000000000-A3JHG_1_1101_3758_10627" 
#>   [71] "M00967_43_000000000-A3JHG_1_1101_4190_14354" 
#>   [72] "M00967_43_000000000-A3JHG_1_1101_4346_14940" 
#>   [73] "M00967_43_000000000-A3JHG_1_1101_5922_6513"  
#>   [74] "M00967_43_000000000-A3JHG_1_1101_6770_21083" 
#>   [75] "M00967_43_000000000-A3JHG_1_1101_6836_23417" 
#>   [76] "M00967_43_000000000-A3JHG_1_1101_6973_14783" 
#>   [77] "M00967_43_000000000-A3JHG_1_1101_7373_14922" 
#>   [78] "M00967_43_000000000-A3JHG_1_1101_7971_11346" 
#>   [79] "M00967_43_000000000-A3JHG_1_1101_8537_7072"  
#>   [80] "M00967_43_000000000-A3JHG_1_1101_8786_10844" 
#>   [81] "M00967_43_000000000-A3JHG_1_1101_8868_4602"  
#>   [82] "M00967_43_000000000-A3JHG_1_1101_8877_3935"  
#>   [83] "M00967_43_000000000-A3JHG_1_1101_8995_24296" 
#>   [84] "M00967_43_000000000-A3JHG_1_1101_9065_10633" 
#>   [85] "M00967_43_000000000-A3JHG_1_1101_9331_5806"  
#>   [86] "M00967_43_000000000-A3JHG_1_1101_9360_18060" 
#>   [87] "M00967_43_000000000-A3JHG_1_1101_9553_14094" 
#>   [88] "M00967_43_000000000-A3JHG_1_1101_9620_19745" 
#>   [89] "M00967_43_000000000-A3JHG_1_1102_10298_5101" 
#>   [90] "M00967_43_000000000-A3JHG_1_1102_10390_20786"
#>   [91] "M00967_43_000000000-A3JHG_1_1102_11115_19075"
#>   [92] "M00967_43_000000000-A3JHG_1_1102_11158_17229"
#>   [93] "M00967_43_000000000-A3JHG_1_1102_11222_12600"
#>   [94] "M00967_43_000000000-A3JHG_1_1102_11324_6735" 
#>   [95] "M00967_43_000000000-A3JHG_1_1102_11664_24324"
#>   [96] "M00967_43_000000000-A3JHG_1_1102_12417_5953" 
#>   [97] "M00967_43_000000000-A3JHG_1_1102_12472_26307"
#>   [98] "M00967_43_000000000-A3JHG_1_1102_12551_11293"
#>   [99] "M00967_43_000000000-A3JHG_1_1102_12689_8055" 
#>  [100] "M00967_43_000000000-A3JHG_1_1102_12944_14923"
#>  [101] "M00967_43_000000000-A3JHG_1_1102_13319_18244"
#>  [102] "M00967_43_000000000-A3JHG_1_1102_13343_4567" 
#>  [103] "M00967_43_000000000-A3JHG_1_1102_14321_10096"
#>  [104] "M00967_43_000000000-A3JHG_1_1102_14586_9284" 
#>  [105] "M00967_43_000000000-A3JHG_1_1102_14615_8833" 
#>  [106] "M00967_43_000000000-A3JHG_1_1102_14900_7576" 
#>  [107] "M00967_43_000000000-A3JHG_1_1102_14905_6542" 
#>  [108] "M00967_43_000000000-A3JHG_1_1102_14928_26203"
#>  [109] "M00967_43_000000000-A3JHG_1_1102_14993_16910"
#>  [110] "M00967_43_000000000-A3JHG_1_1102_15058_11924"
#>  [111] "M00967_43_000000000-A3JHG_1_1102_15350_5623" 
#>  [112] "M00967_43_000000000-A3JHG_1_1102_15600_19124"
#>  [113] "M00967_43_000000000-A3JHG_1_1102_15632_13472"
#>  [114] "M00967_43_000000000-A3JHG_1_1102_16785_23687"
#>  [115] "M00967_43_000000000-A3JHG_1_1102_16964_23052"
#>  [116] "M00967_43_000000000-A3JHG_1_1102_17763_13764"
#>  [117] "M00967_43_000000000-A3JHG_1_1102_17902_26485"
#>  [118] "M00967_43_000000000-A3JHG_1_1102_18244_26058"
#>  [119] "M00967_43_000000000-A3JHG_1_1102_18640_14309"
#>  [120] "M00967_43_000000000-A3JHG_1_1102_18759_22050"
#>  [121] "M00967_43_000000000-A3JHG_1_1102_19038_11166"
#>  [122] "M00967_43_000000000-A3JHG_1_1102_19130_22074"
#>  [123] "M00967_43_000000000-A3JHG_1_1102_19172_4834" 
#>  [124] "M00967_43_000000000-A3JHG_1_1102_19764_24086"
#>  [125] "M00967_43_000000000-A3JHG_1_1102_20639_13713"
#>  [126] "M00967_43_000000000-A3JHG_1_1102_20738_4913" 
#>  [127] "M00967_43_000000000-A3JHG_1_1102_20961_15950"
#>  [128] "M00967_43_000000000-A3JHG_1_1102_2114_15227" 
#>  [129] "M00967_43_000000000-A3JHG_1_1102_21170_25493"
#>  [130] "M00967_43_000000000-A3JHG_1_1102_21598_16589"
#>  [131] "M00967_43_000000000-A3JHG_1_1102_21808_10504"
#>  [132] "M00967_43_000000000-A3JHG_1_1102_22477_7326" 
#>  [133] "M00967_43_000000000-A3JHG_1_1102_22608_17863"
#>  [134] "M00967_43_000000000-A3JHG_1_1102_22955_21067"
#>  [135] "M00967_43_000000000-A3JHG_1_1102_24758_6643" 
#>  [136] "M00967_43_000000000-A3JHG_1_1102_24841_18959"
#>  [137] "M00967_43_000000000-A3JHG_1_1102_25331_23286"
#>  [138] "M00967_43_000000000-A3JHG_1_1102_26098_19093"
#>  [139] "M00967_43_000000000-A3JHG_1_1102_26308_14496"
#>  [140] "M00967_43_000000000-A3JHG_1_1102_26403_8487" 
#>  [141] "M00967_43_000000000-A3JHG_1_1102_26826_15050"
#>  [142] "M00967_43_000000000-A3JHG_1_1102_27180_15478"
#>  [143] "M00967_43_000000000-A3JHG_1_1102_27414_13728"
#>  [144] "M00967_43_000000000-A3JHG_1_1102_27837_16697"
#>  [145] "M00967_43_000000000-A3JHG_1_1102_27953_19313"
#>  [146] "M00967_43_000000000-A3JHG_1_1102_3841_17770" 
#>  [147] "M00967_43_000000000-A3JHG_1_1102_3890_17221" 
#>  [148] "M00967_43_000000000-A3JHG_1_1102_4170_17056" 
#>  [149] "M00967_43_000000000-A3JHG_1_1102_4962_7612"  
#>  [150] "M00967_43_000000000-A3JHG_1_1102_5327_17123" 
#>  [151] "M00967_43_000000000-A3JHG_1_1102_5497_9039"  
#>  [152] "M00967_43_000000000-A3JHG_1_1102_5575_6810"  
#>  [153] "M00967_43_000000000-A3JHG_1_1102_6774_6343"  
#>  [154] "M00967_43_000000000-A3JHG_1_1102_7684_6622"  
#>  [155] "M00967_43_000000000-A3JHG_1_1102_7847_17999" 
#>  [156] "M00967_43_000000000-A3JHG_1_1102_8406_20325" 
#>  [157] "M00967_43_000000000-A3JHG_1_1102_8524_16122" 
#>  [158] "M00967_43_000000000-A3JHG_1_1102_8796_20607" 
#>  [159] "M00967_43_000000000-A3JHG_1_1102_8887_21995" 
#>  [160] "M00967_43_000000000-A3JHG_1_1102_9244_9305"  
#>  [161] "M00967_43_000000000-A3JHG_1_1103_10339_19353"
#>  [162] "M00967_43_000000000-A3JHG_1_1103_10830_22855"
#>  [163] "M00967_43_000000000-A3JHG_1_1103_10837_19652"
#>  [164] "M00967_43_000000000-A3JHG_1_1103_11367_13197"
#>  [165] "M00967_43_000000000-A3JHG_1_1103_12046_28270"
#>  [166] "M00967_43_000000000-A3JHG_1_1103_12195_22317"
#>  [167] "M00967_43_000000000-A3JHG_1_1103_12270_15236"
#>  [168] "M00967_43_000000000-A3JHG_1_1103_12442_25097"
#>  [169] "M00967_43_000000000-A3JHG_1_1103_12856_12076"
#>  [170] "M00967_43_000000000-A3JHG_1_1103_12955_24887"
#>  [171] "M00967_43_000000000-A3JHG_1_1103_13035_5452" 
#>  [172] "M00967_43_000000000-A3JHG_1_1103_13340_15811"
#>  [173] "M00967_43_000000000-A3JHG_1_1103_13364_5496" 
#>  [174] "M00967_43_000000000-A3JHG_1_1103_13966_3813" 
#>  [175] "M00967_43_000000000-A3JHG_1_1103_14435_10999"
#>  [176] "M00967_43_000000000-A3JHG_1_1103_14518_16099"
#>  [177] "M00967_43_000000000-A3JHG_1_1103_14890_4136" 
#>  [178] "M00967_43_000000000-A3JHG_1_1103_15776_7109" 
#>  [179] "M00967_43_000000000-A3JHG_1_1103_16201_28240"
#>  [180] "M00967_43_000000000-A3JHG_1_1103_16766_20121"
#>  [181] "M00967_43_000000000-A3JHG_1_1103_18369_27790"
#>  [182] "M00967_43_000000000-A3JHG_1_1103_18432_8178" 
#>  [183] "M00967_43_000000000-A3JHG_1_1103_18945_9651" 
#>  [184] "M00967_43_000000000-A3JHG_1_1103_19321_18100"
#>  [185] "M00967_43_000000000-A3JHG_1_1103_19411_10296"
#>  [186] "M00967_43_000000000-A3JHG_1_1103_19447_15932"
#>  [187] "M00967_43_000000000-A3JHG_1_1103_19870_21567"
#>  [188] "M00967_43_000000000-A3JHG_1_1103_20010_17589"
#>  [189] "M00967_43_000000000-A3JHG_1_1103_21051_5371" 
#>  [190] "M00967_43_000000000-A3JHG_1_1103_22446_9538" 
#>  [191] "M00967_43_000000000-A3JHG_1_1103_22490_21890"
#>  [192] "M00967_43_000000000-A3JHG_1_1103_22605_27014"
#>  [193] "M00967_43_000000000-A3JHG_1_1103_22654_26891"
#>  [194] "M00967_43_000000000-A3JHG_1_1103_22674_6370" 
#>  [195] "M00967_43_000000000-A3JHG_1_1103_22686_12139"
#>  [196] "M00967_43_000000000-A3JHG_1_1103_22791_24807"
#>  [197] "M00967_43_000000000-A3JHG_1_1103_22842_22813"
#>  [198] "M00967_43_000000000-A3JHG_1_1103_22959_4225" 
#>  [199] "M00967_43_000000000-A3JHG_1_1103_23102_6138" 
#>  [200] "M00967_43_000000000-A3JHG_1_1103_23279_25937"
#>  [201] "M00967_43_000000000-A3JHG_1_1103_2359_11998" 
#>  [202] "M00967_43_000000000-A3JHG_1_1103_23742_6714" 
#>  [203] "M00967_43_000000000-A3JHG_1_1103_23908_15902"
#>  [204] "M00967_43_000000000-A3JHG_1_1103_2408_17510" 
#>  [205] "M00967_43_000000000-A3JHG_1_1103_24229_16660"
#>  [206] "M00967_43_000000000-A3JHG_1_1103_24387_17630"
#>  [207] "M00967_43_000000000-A3JHG_1_1103_24429_8227" 
#>  [208] "M00967_43_000000000-A3JHG_1_1103_24460_16041"
#>  [209] "M00967_43_000000000-A3JHG_1_1103_25888_23585"
#>  [210] "M00967_43_000000000-A3JHG_1_1103_26302_22051"
#>  [211] "M00967_43_000000000-A3JHG_1_1103_26542_17589"
#>  [212] "M00967_43_000000000-A3JHG_1_1103_26580_14708"
#>  [213] "M00967_43_000000000-A3JHG_1_1103_27900_10832"
#>  [214] "M00967_43_000000000-A3JHG_1_1103_28341_20456"
#>  [215] "M00967_43_000000000-A3JHG_1_1103_28451_19405"
#>  [216] "M00967_43_000000000-A3JHG_1_1103_3710_21770" 
#>  [217] "M00967_43_000000000-A3JHG_1_1103_3864_17599" 
#>  [218] "M00967_43_000000000-A3JHG_1_1103_4808_12955" 
#>  [219] "M00967_43_000000000-A3JHG_1_1103_5171_14027" 
#>  [220] "M00967_43_000000000-A3JHG_1_1103_5501_16588" 
#>  [221] "M00967_43_000000000-A3JHG_1_1103_5754_14689" 
#>  [222] "M00967_43_000000000-A3JHG_1_1103_6362_13340" 
#>  [223] "M00967_43_000000000-A3JHG_1_1103_6514_8778"  
#>  [224] "M00967_43_000000000-A3JHG_1_1103_6808_22767" 
#>  [225] "M00967_43_000000000-A3JHG_1_1103_7049_5487"  
#>  [226] "M00967_43_000000000-A3JHG_1_1103_7754_25974" 
#>  [227] "M00967_43_000000000-A3JHG_1_1103_7866_21570" 
#>  [228] "M00967_43_000000000-A3JHG_1_1103_7913_12679" 
#>  [229] "M00967_43_000000000-A3JHG_1_1103_7927_15626" 
#>  [230] "M00967_43_000000000-A3JHG_1_1103_8349_14557" 
#>  [231] "M00967_43_000000000-A3JHG_1_1103_8361_9953"  
#>  [232] "M00967_43_000000000-A3JHG_1_1103_8412_8356"  
#>  [233] "M00967_43_000000000-A3JHG_1_1103_8900_10891" 
#>  [234] "M00967_43_000000000-A3JHG_1_1103_9123_12021" 
#>  [235] "M00967_43_000000000-A3JHG_1_1103_9133_12955" 
#>  [236] "M00967_43_000000000-A3JHG_1_1103_9570_21383" 
#>  [237] "M00967_43_000000000-A3JHG_1_1103_9802_9756"  
#>  [238] "M00967_43_000000000-A3JHG_1_1103_9913_21564" 
#>  [239] "M00967_43_000000000-A3JHG_1_1104_10411_22472"
#>  [240] "M00967_43_000000000-A3JHG_1_1104_10691_10850"
#>  [241] "M00967_43_000000000-A3JHG_1_1104_10756_24942"
#>  [242] "M00967_43_000000000-A3JHG_1_1104_10834_8009" 
#>  [243] "M00967_43_000000000-A3JHG_1_1104_11335_21750"
#>  [244] "M00967_43_000000000-A3JHG_1_1104_11535_9221" 
#>  [245] "M00967_43_000000000-A3JHG_1_1104_11751_17981"
#>  [246] "M00967_43_000000000-A3JHG_1_1104_11900_6285" 
#>  [247] "M00967_43_000000000-A3JHG_1_1104_12015_3840" 
#>  [248] "M00967_43_000000000-A3JHG_1_1104_12197_8664" 
#>  [249] "M00967_43_000000000-A3JHG_1_1104_12222_5327" 
#>  [250] "M00967_43_000000000-A3JHG_1_1104_12675_27249"
#>  [251] "M00967_43_000000000-A3JHG_1_1104_13037_4245" 
#>  [252] "M00967_43_000000000-A3JHG_1_1104_13663_10565"
#>  [253] "M00967_43_000000000-A3JHG_1_1104_13776_17053"
#>  [254] "M00967_43_000000000-A3JHG_1_1104_14213_12414"
#>  [255] "M00967_43_000000000-A3JHG_1_1104_14252_27883"
#>  [256] "M00967_43_000000000-A3JHG_1_1104_14415_16643"
#>  [257] "M00967_43_000000000-A3JHG_1_1104_15033_26207"
#>  [258] "M00967_43_000000000-A3JHG_1_1104_15283_26901"
#>  [259] "M00967_43_000000000-A3JHG_1_1104_15705_5052" 
#>  [260] "M00967_43_000000000-A3JHG_1_1104_15748_9889" 
#>  [261] "M00967_43_000000000-A3JHG_1_1104_15789_5235" 
#>  [262] "M00967_43_000000000-A3JHG_1_1104_17110_15247"
#>  [263] "M00967_43_000000000-A3JHG_1_1104_17143_19056"
#>  [264] "M00967_43_000000000-A3JHG_1_1104_17151_24507"
#>  [265] "M00967_43_000000000-A3JHG_1_1104_17583_11048"
#>  [266] "M00967_43_000000000-A3JHG_1_1104_17756_12554"
#>  [267] "M00967_43_000000000-A3JHG_1_1104_18175_14208"
#>  [268] "M00967_43_000000000-A3JHG_1_1104_18319_13388"
#>  [269] "M00967_43_000000000-A3JHG_1_1104_18924_18211"
#>  [270] "M00967_43_000000000-A3JHG_1_1104_19278_14733"
#>  [271] "M00967_43_000000000-A3JHG_1_1104_19441_6353" 
#>  [272] "M00967_43_000000000-A3JHG_1_1104_19580_19072"
#>  [273] "M00967_43_000000000-A3JHG_1_1104_19589_20033"
#>  [274] "M00967_43_000000000-A3JHG_1_1104_19763_16198"
#>  [275] "M00967_43_000000000-A3JHG_1_1104_20104_23610"
#>  [276] "M00967_43_000000000-A3JHG_1_1104_20303_14415"
#>  [277] "M00967_43_000000000-A3JHG_1_1104_20592_3815" 
#>  [278] "M00967_43_000000000-A3JHG_1_1104_20673_8839" 
#>  [279] "M00967_43_000000000-A3JHG_1_1104_20877_7748" 
#>  [280] "M00967_43_000000000-A3JHG_1_1104_21881_13416"
#>  [281] "M00967_43_000000000-A3JHG_1_1104_22299_24998"
#>  [282] "M00967_43_000000000-A3JHG_1_1104_23236_5434" 
#>  [283] "M00967_43_000000000-A3JHG_1_1104_23248_7732" 
#>  [284] "M00967_43_000000000-A3JHG_1_1104_23689_6395" 
#>  [285] "M00967_43_000000000-A3JHG_1_1104_23719_15777"
#>  [286] "M00967_43_000000000-A3JHG_1_1104_23849_8759" 
#>  [287] "M00967_43_000000000-A3JHG_1_1104_23895_17769"
#>  [288] "M00967_43_000000000-A3JHG_1_1104_24088_7145" 
#>  [289] "M00967_43_000000000-A3JHG_1_1104_24119_21632"
#>  [290] "M00967_43_000000000-A3JHG_1_1104_24270_11650"
#>  [291] "M00967_43_000000000-A3JHG_1_1104_24923_8045" 
#>  [292] "M00967_43_000000000-A3JHG_1_1104_24989_12614"
#>  [293] "M00967_43_000000000-A3JHG_1_1104_25277_11666"
#>  [294] "M00967_43_000000000-A3JHG_1_1104_25962_6708" 
#>  [295] "M00967_43_000000000-A3JHG_1_1104_27270_10834"
#>  [296] "M00967_43_000000000-A3JHG_1_1104_27578_17650"
#>  [297] "M00967_43_000000000-A3JHG_1_1104_27702_14915"
#>  [298] "M00967_43_000000000-A3JHG_1_1104_29281_13264"
#>  [299] "M00967_43_000000000-A3JHG_1_1104_4327_10073" 
#>  [300] "M00967_43_000000000-A3JHG_1_1104_5274_12644" 
#>  [301] "M00967_43_000000000-A3JHG_1_1104_5458_8163"  
#>  [302] "M00967_43_000000000-A3JHG_1_1104_5923_11043" 
#>  [303] "M00967_43_000000000-A3JHG_1_1104_6527_16493" 
#>  [304] "M00967_43_000000000-A3JHG_1_1104_7070_12571" 
#>  [305] "M00967_43_000000000-A3JHG_1_1104_7366_7146"  
#>  [306] "M00967_43_000000000-A3JHG_1_1104_7613_8427"  
#>  [307] "M00967_43_000000000-A3JHG_1_1104_7782_10984" 
#>  [308] "M00967_43_000000000-A3JHG_1_1104_8057_17803" 
#>  [309] "M00967_43_000000000-A3JHG_1_1104_8585_16521" 
#>  [310] "M00967_43_000000000-A3JHG_1_1104_9048_17067" 
#>  [311] "M00967_43_000000000-A3JHG_1_1105_10343_20641"
#>  [312] "M00967_43_000000000-A3JHG_1_1105_10652_23248"
#>  [313] "M00967_43_000000000-A3JHG_1_1105_11157_16891"
#>  [314] "M00967_43_000000000-A3JHG_1_1105_11242_27278"
#>  [315] "M00967_43_000000000-A3JHG_1_1105_11497_20872"
#>  [316] "M00967_43_000000000-A3JHG_1_1105_11964_4686" 
#>  [317] "M00967_43_000000000-A3JHG_1_1105_11983_14742"
#>  [318] "M00967_43_000000000-A3JHG_1_1105_12311_16407"
#>  [319] "M00967_43_000000000-A3JHG_1_1105_12460_19617"
#>  [320] "M00967_43_000000000-A3JHG_1_1105_12761_22353"
#>  [321] "M00967_43_000000000-A3JHG_1_1105_13025_22602"
#>  [322] "M00967_43_000000000-A3JHG_1_1105_13657_20414"
#>  [323] "M00967_43_000000000-A3JHG_1_1105_14169_25963"
#>  [324] "M00967_43_000000000-A3JHG_1_1105_14547_12843"
#>  [325] "M00967_43_000000000-A3JHG_1_1105_14835_21938"
#>  [326] "M00967_43_000000000-A3JHG_1_1105_14922_22316"
#>  [327] "M00967_43_000000000-A3JHG_1_1105_14935_17739"
#>  [328] "M00967_43_000000000-A3JHG_1_1105_14971_7732" 
#>  [329] "M00967_43_000000000-A3JHG_1_1105_15090_18340"
#>  [330] "M00967_43_000000000-A3JHG_1_1105_15318_9427" 
#>  [331] "M00967_43_000000000-A3JHG_1_1105_15500_1801" 
#>  [332] "M00967_43_000000000-A3JHG_1_1105_15558_26326"
#>  [333] "M00967_43_000000000-A3JHG_1_1105_15601_18021"
#>  [334] "M00967_43_000000000-A3JHG_1_1105_15636_27192"
#>  [335] "M00967_43_000000000-A3JHG_1_1105_15675_15142"
#>  [336] "M00967_43_000000000-A3JHG_1_1105_15878_26239"
#>  [337] "M00967_43_000000000-A3JHG_1_1105_16167_11796"
#>  [338] "M00967_43_000000000-A3JHG_1_1105_16224_23592"
#>  [339] "M00967_43_000000000-A3JHG_1_1105_16333_26309"
#>  [340] "M00967_43_000000000-A3JHG_1_1105_16624_8157" 
#>  [341] "M00967_43_000000000-A3JHG_1_1105_17600_8589" 
#>  [342] "M00967_43_000000000-A3JHG_1_1105_17602_8360" 
#>  [343] "M00967_43_000000000-A3JHG_1_1105_17706_11167"
#>  [344] "M00967_43_000000000-A3JHG_1_1105_18275_26860"
#>  [345] "M00967_43_000000000-A3JHG_1_1105_18475_22342"
#>  [346] "M00967_43_000000000-A3JHG_1_1105_18566_14540"
#>  [347] "M00967_43_000000000-A3JHG_1_1105_18918_5077" 
#>  [348] "M00967_43_000000000-A3JHG_1_1105_18988_15558"
#>  [349] "M00967_43_000000000-A3JHG_1_1105_19059_19452"
#>  [350] "M00967_43_000000000-A3JHG_1_1105_19259_21539"
#>  [351] "M00967_43_000000000-A3JHG_1_1105_20176_26018"
#>  [352] "M00967_43_000000000-A3JHG_1_1105_20464_8555" 
#>  [353] "M00967_43_000000000-A3JHG_1_1105_20465_8611" 
#>  [354] "M00967_43_000000000-A3JHG_1_1105_20476_12712"
#>  [355] "M00967_43_000000000-A3JHG_1_1105_21018_4980" 
#>  [356] "M00967_43_000000000-A3JHG_1_1105_21138_11762"
#>  [357] "M00967_43_000000000-A3JHG_1_1105_22037_19357"
#>  [358] "M00967_43_000000000-A3JHG_1_1105_22101_17039"
#>  [359] "M00967_43_000000000-A3JHG_1_1105_22508_16792"
#>  [360] "M00967_43_000000000-A3JHG_1_1105_22535_20381"
#>  [361] "M00967_43_000000000-A3JHG_1_1105_22581_18347"
#>  [362] "M00967_43_000000000-A3JHG_1_1105_23025_12058"
#>  [363] "M00967_43_000000000-A3JHG_1_1105_23025_12099"
#>  [364] "M00967_43_000000000-A3JHG_1_1105_23153_7273" 
#>  [365] "M00967_43_000000000-A3JHG_1_1105_23231_4481" 
#>  [366] "M00967_43_000000000-A3JHG_1_1105_23612_18072"
#>  [367] "M00967_43_000000000-A3JHG_1_1105_23740_10259"
#>  [368] "M00967_43_000000000-A3JHG_1_1105_23929_25607"
#>  [369] "M00967_43_000000000-A3JHG_1_1105_24853_19323"
#>  [370] "M00967_43_000000000-A3JHG_1_1105_25210_11850"
#>  [371] "M00967_43_000000000-A3JHG_1_1105_25410_11631"
#>  [372] "M00967_43_000000000-A3JHG_1_1105_25421_15745"
#>  [373] "M00967_43_000000000-A3JHG_1_1105_25554_16982"
#>  [374] "M00967_43_000000000-A3JHG_1_1105_25642_17588"
#>  [375] "M00967_43_000000000-A3JHG_1_1105_25694_18609"
#>  [376] "M00967_43_000000000-A3JHG_1_1105_26683_22024"
#>  [377] "M00967_43_000000000-A3JHG_1_1105_27296_22083"
#>  [378] "M00967_43_000000000-A3JHG_1_1105_27543_11258"
#>  [379] "M00967_43_000000000-A3JHG_1_1105_27806_20120"
#>  [380] "M00967_43_000000000-A3JHG_1_1105_3088_15065" 
#>  [381] "M00967_43_000000000-A3JHG_1_1105_5158_15329" 
#>  [382] "M00967_43_000000000-A3JHG_1_1105_5457_15684" 
#>  [383] "M00967_43_000000000-A3JHG_1_1105_5882_22900" 
#>  [384] "M00967_43_000000000-A3JHG_1_1105_6908_19395" 
#>  [385] "M00967_43_000000000-A3JHG_1_1105_7287_9945"  
#>  [386] "M00967_43_000000000-A3JHG_1_1105_7593_10522" 
#>  [387] "M00967_43_000000000-A3JHG_1_1105_8356_13589" 
#>  [388] "M00967_43_000000000-A3JHG_1_1105_8707_15634" 
#>  [389] "M00967_43_000000000-A3JHG_1_1105_8938_25075" 
#>  [390] "M00967_43_000000000-A3JHG_1_1105_9050_25130" 
#>  [391] "M00967_43_000000000-A3JHG_1_1105_9179_20196" 
#>  [392] "M00967_43_000000000-A3JHG_1_1105_9374_22772" 
#>  [393] "M00967_43_000000000-A3JHG_1_1105_9515_26311" 
#>  [394] "M00967_43_000000000-A3JHG_1_1105_9868_4758"  
#>  [395] "M00967_43_000000000-A3JHG_1_1105_9992_20296" 
#>  [396] "M00967_43_000000000-A3JHG_1_1106_10456_18286"
#>  [397] "M00967_43_000000000-A3JHG_1_1106_10890_9512" 
#>  [398] "M00967_43_000000000-A3JHG_1_1106_10941_17222"
#>  [399] "M00967_43_000000000-A3JHG_1_1106_11032_21384"
#>  [400] "M00967_43_000000000-A3JHG_1_1106_11208_20468"
#>  [401] "M00967_43_000000000-A3JHG_1_1106_11216_2465" 
#>  [402] "M00967_43_000000000-A3JHG_1_1106_11240_22282"
#>  [403] "M00967_43_000000000-A3JHG_1_1106_11283_15681"
#>  [404] "M00967_43_000000000-A3JHG_1_1106_11865_6449" 
#>  [405] "M00967_43_000000000-A3JHG_1_1106_12231_13452"
#>  [406] "M00967_43_000000000-A3JHG_1_1106_12236_14361"
#>  [407] "M00967_43_000000000-A3JHG_1_1106_12867_4666" 
#>  [408] "M00967_43_000000000-A3JHG_1_1106_14252_8640" 
#>  [409] "M00967_43_000000000-A3JHG_1_1106_14254_12461"
#>  [410] "M00967_43_000000000-A3JHG_1_1106_14256_14203"
#>  [411] "M00967_43_000000000-A3JHG_1_1106_14915_20467"
#>  [412] "M00967_43_000000000-A3JHG_1_1106_14993_9292" 
#>  [413] "M00967_43_000000000-A3JHG_1_1106_15774_18346"
#>  [414] "M00967_43_000000000-A3JHG_1_1106_15894_8108" 
#>  [415] "M00967_43_000000000-A3JHG_1_1106_15955_6621" 
#>  [416] "M00967_43_000000000-A3JHG_1_1106_15986_21920"
#>  [417] "M00967_43_000000000-A3JHG_1_1106_16049_12791"
#>  [418] "M00967_43_000000000-A3JHG_1_1106_16075_4721" 
#>  [419] "M00967_43_000000000-A3JHG_1_1106_17225_7144" 
#>  [420] "M00967_43_000000000-A3JHG_1_1106_17565_8490" 
#>  [421] "M00967_43_000000000-A3JHG_1_1106_17680_15930"
#>  [422] "M00967_43_000000000-A3JHG_1_1106_17692_2367" 
#>  [423] "M00967_43_000000000-A3JHG_1_1106_17965_25181"
#>  [424] "M00967_43_000000000-A3JHG_1_1106_18033_2621" 
#>  [425] "M00967_43_000000000-A3JHG_1_1106_18148_6365" 
#>  [426] "M00967_43_000000000-A3JHG_1_1106_18589_17554"
#>  [427] "M00967_43_000000000-A3JHG_1_1106_18846_19788"
#>  [428] "M00967_43_000000000-A3JHG_1_1106_19006_23053"
#>  [429] "M00967_43_000000000-A3JHG_1_1106_19617_4733" 
#>  [430] "M00967_43_000000000-A3JHG_1_1106_19777_24988"
#>  [431] "M00967_43_000000000-A3JHG_1_1106_20020_16264"
#>  [432] "M00967_43_000000000-A3JHG_1_1106_20132_15845"
#>  [433] "M00967_43_000000000-A3JHG_1_1106_20287_11044"
#>  [434] "M00967_43_000000000-A3JHG_1_1106_20644_25868"
#>  [435] "M00967_43_000000000-A3JHG_1_1106_20737_14683"
#>  [436] "M00967_43_000000000-A3JHG_1_1106_20923_14281"
#>  [437] "M00967_43_000000000-A3JHG_1_1106_21375_27109"
#>  [438] "M00967_43_000000000-A3JHG_1_1106_21486_19228"
#>  [439] "M00967_43_000000000-A3JHG_1_1106_22389_22397"
#>  [440] "M00967_43_000000000-A3JHG_1_1106_22705_6123" 
#>  [441] "M00967_43_000000000-A3JHG_1_1106_22826_4799" 
#>  [442] "M00967_43_000000000-A3JHG_1_1106_22881_7086" 
#>  [443] "M00967_43_000000000-A3JHG_1_1106_2293_12053" 
#>  [444] "M00967_43_000000000-A3JHG_1_1106_23059_17918"
#>  [445] "M00967_43_000000000-A3JHG_1_1106_23133_8169" 
#>  [446] "M00967_43_000000000-A3JHG_1_1106_23262_25144"
#>  [447] "M00967_43_000000000-A3JHG_1_1106_23455_4853" 
#>  [448] "M00967_43_000000000-A3JHG_1_1106_23502_7683" 
#>  [449] "M00967_43_000000000-A3JHG_1_1106_23697_7097" 
#>  [450] "M00967_43_000000000-A3JHG_1_1106_23989_25110"
#>  [451] "M00967_43_000000000-A3JHG_1_1106_24036_10250"
#>  [452] "M00967_43_000000000-A3JHG_1_1106_24084_10563"
#>  [453] "M00967_43_000000000-A3JHG_1_1106_24343_18277"
#>  [454] "M00967_43_000000000-A3JHG_1_1106_24378_9680" 
#>  [455] "M00967_43_000000000-A3JHG_1_1106_24555_20879"
#>  [456] "M00967_43_000000000-A3JHG_1_1106_25034_19107"
#>  [457] "M00967_43_000000000-A3JHG_1_1106_25330_8604" 
#>  [458] "M00967_43_000000000-A3JHG_1_1106_25389_22814"
#>  [459] "M00967_43_000000000-A3JHG_1_1106_25828_15956"
#>  [460] "M00967_43_000000000-A3JHG_1_1106_25979_21288"
#>  [461] "M00967_43_000000000-A3JHG_1_1106_26087_11342"
#>  [462] "M00967_43_000000000-A3JHG_1_1106_26352_22378"
#>  [463] "M00967_43_000000000-A3JHG_1_1106_26592_14233"
#>  [464] "M00967_43_000000000-A3JHG_1_1106_26888_10548"
#>  [465] "M00967_43_000000000-A3JHG_1_1106_27997_12482"
#>  [466] "M00967_43_000000000-A3JHG_1_1106_28442_12734"
#>  [467] "M00967_43_000000000-A3JHG_1_1106_28698_12913"
#>  [468] "M00967_43_000000000-A3JHG_1_1106_3789_8994"  
#>  [469] "M00967_43_000000000-A3JHG_1_1106_4308_15219" 
#>  [470] "M00967_43_000000000-A3JHG_1_1106_5336_9821"  
#>  [471] "M00967_43_000000000-A3JHG_1_1106_5459_18973" 
#>  [472] "M00967_43_000000000-A3JHG_1_1106_5527_16685" 
#>  [473] "M00967_43_000000000-A3JHG_1_1106_5620_18651" 
#>  [474] "M00967_43_000000000-A3JHG_1_1106_6004_8081"  
#>  [475] "M00967_43_000000000-A3JHG_1_1106_6775_23397" 
#>  [476] "M00967_43_000000000-A3JHG_1_1106_7344_20185" 
#>  [477] "M00967_43_000000000-A3JHG_1_1106_7651_17094" 
#>  [478] "M00967_43_000000000-A3JHG_1_1106_7807_21999" 
#>  [479] "M00967_43_000000000-A3JHG_1_1106_8066_16367" 
#>  [480] "M00967_43_000000000-A3JHG_1_1106_8166_23839" 
#>  [481] "M00967_43_000000000-A3JHG_1_1106_8183_16295" 
#>  [482] "M00967_43_000000000-A3JHG_1_1106_8303_10880" 
#>  [483] "M00967_43_000000000-A3JHG_1_1106_8681_8228"  
#>  [484] "M00967_43_000000000-A3JHG_1_1106_8718_15420" 
#>  [485] "M00967_43_000000000-A3JHG_1_1106_9147_22015" 
#>  [486] "M00967_43_000000000-A3JHG_1_1106_9232_9214"  
#>  [487] "M00967_43_000000000-A3JHG_1_1106_9558_4409"  
#>  [488] "M00967_43_000000000-A3JHG_1_1106_9733_12197" 
#>  [489] "M00967_43_000000000-A3JHG_1_1107_10545_13535"
#>  [490] "M00967_43_000000000-A3JHG_1_1107_10632_26179"
#>  [491] "M00967_43_000000000-A3JHG_1_1107_10661_18652"
#>  [492] "M00967_43_000000000-A3JHG_1_1107_10811_19778"
#>  [493] "M00967_43_000000000-A3JHG_1_1107_11459_12453"
#>  [494] "M00967_43_000000000-A3JHG_1_1107_11517_8173" 
#>  [495] "M00967_43_000000000-A3JHG_1_1107_12586_26826"
#>  [496] "M00967_43_000000000-A3JHG_1_1107_12904_20713"
#>  [497] "M00967_43_000000000-A3JHG_1_1107_13057_19114"
#>  [498] "M00967_43_000000000-A3JHG_1_1107_13116_15821"
#>  [499] "M00967_43_000000000-A3JHG_1_1107_13192_19802"
#>  [500] "M00967_43_000000000-A3JHG_1_1107_14088_12072"
#>  [501] "M00967_43_000000000-A3JHG_1_1107_14329_23174"
#>  [502] "M00967_43_000000000-A3JHG_1_1107_14437_24723"
#>  [503] "M00967_43_000000000-A3JHG_1_1107_14488_13179"
#>  [504] "M00967_43_000000000-A3JHG_1_1107_14548_18902"
#>  [505] "M00967_43_000000000-A3JHG_1_1107_14868_16744"
#>  [506] "M00967_43_000000000-A3JHG_1_1107_15327_17780"
#>  [507] "M00967_43_000000000-A3JHG_1_1107_15441_15657"
#>  [508] "M00967_43_000000000-A3JHG_1_1107_15635_4811" 
#>  [509] "M00967_43_000000000-A3JHG_1_1107_15708_17981"
#>  [510] "M00967_43_000000000-A3JHG_1_1107_15750_18592"
#>  [511] "M00967_43_000000000-A3JHG_1_1107_15765_25536"
#>  [512] "M00967_43_000000000-A3JHG_1_1107_16297_12364"
#>  [513] "M00967_43_000000000-A3JHG_1_1107_16665_12409"
#>  [514] "M00967_43_000000000-A3JHG_1_1107_17032_8072" 
#>  [515] "M00967_43_000000000-A3JHG_1_1107_17052_3882" 
#>  [516] "M00967_43_000000000-A3JHG_1_1107_17173_6497" 
#>  [517] "M00967_43_000000000-A3JHG_1_1107_18594_19714"
#>  [518] "M00967_43_000000000-A3JHG_1_1107_18662_2044" 
#>  [519] "M00967_43_000000000-A3JHG_1_1107_18879_6171" 
#>  [520] "M00967_43_000000000-A3JHG_1_1107_18906_15466"
#>  [521] "M00967_43_000000000-A3JHG_1_1107_19005_21895"
#>  [522] "M00967_43_000000000-A3JHG_1_1107_19088_18601"
#>  [523] "M00967_43_000000000-A3JHG_1_1107_19357_22811"
#>  [524] "M00967_43_000000000-A3JHG_1_1107_19562_20695"
#>  [525] "M00967_43_000000000-A3JHG_1_1107_20108_5368" 
#>  [526] "M00967_43_000000000-A3JHG_1_1107_20748_16898"
#>  [527] "M00967_43_000000000-A3JHG_1_1107_21038_8569" 
#>  [528] "M00967_43_000000000-A3JHG_1_1107_21132_10931"
#>  [529] "M00967_43_000000000-A3JHG_1_1107_21227_17435"
#>  [530] "M00967_43_000000000-A3JHG_1_1107_21271_3203" 
#>  [531] "M00967_43_000000000-A3JHG_1_1107_21577_23733"
#>  [532] "M00967_43_000000000-A3JHG_1_1107_21734_25420"
#>  [533] "M00967_43_000000000-A3JHG_1_1107_21848_23651"
#>  [534] "M00967_43_000000000-A3JHG_1_1107_22214_8324" 
#>  [535] "M00967_43_000000000-A3JHG_1_1107_22280_22887"
#>  [536] "M00967_43_000000000-A3JHG_1_1107_22428_23154"
#>  [537] "M00967_43_000000000-A3JHG_1_1107_22525_22794"
#>  [538] "M00967_43_000000000-A3JHG_1_1107_22580_21773"
#>  [539] "M00967_43_000000000-A3JHG_1_1107_22659_17429"
#>  [540] "M00967_43_000000000-A3JHG_1_1107_23336_26415"
#>  [541] "M00967_43_000000000-A3JHG_1_1107_23381_22618"
#>  [542] "M00967_43_000000000-A3JHG_1_1107_23638_4659" 
#>  [543] "M00967_43_000000000-A3JHG_1_1107_23853_7456" 
#>  [544] "M00967_43_000000000-A3JHG_1_1107_24649_13577"
#>  [545] "M00967_43_000000000-A3JHG_1_1107_2468_14102" 
#>  [546] "M00967_43_000000000-A3JHG_1_1107_24854_21229"
#>  [547] "M00967_43_000000000-A3JHG_1_1107_24891_17279"
#>  [548] "M00967_43_000000000-A3JHG_1_1107_24903_5978" 
#>  [549] "M00967_43_000000000-A3JHG_1_1107_25867_16609"
#>  [550] "M00967_43_000000000-A3JHG_1_1107_26400_19526"
#>  [551] "M00967_43_000000000-A3JHG_1_1107_26583_8053" 
#>  [552] "M00967_43_000000000-A3JHG_1_1107_26777_15267"
#>  [553] "M00967_43_000000000-A3JHG_1_1107_26939_12098"
#>  [554] "M00967_43_000000000-A3JHG_1_1107_27189_20696"
#>  [555] "M00967_43_000000000-A3JHG_1_1107_27199_15268"
#>  [556] "M00967_43_000000000-A3JHG_1_1107_27450_19498"
#>  [557] "M00967_43_000000000-A3JHG_1_1107_27878_20060"
#>  [558] "M00967_43_000000000-A3JHG_1_1107_3548_11064" 
#>  [559] "M00967_43_000000000-A3JHG_1_1107_4251_18782" 
#>  [560] "M00967_43_000000000-A3JHG_1_1107_4609_18391" 
#>  [561] "M00967_43_000000000-A3JHG_1_1107_5337_17921" 
#>  [562] "M00967_43_000000000-A3JHG_1_1107_6365_24079" 
#>  [563] "M00967_43_000000000-A3JHG_1_1107_6417_21653" 
#>  [564] "M00967_43_000000000-A3JHG_1_1107_6465_14196" 
#>  [565] "M00967_43_000000000-A3JHG_1_1107_7191_15338" 
#>  [566] "M00967_43_000000000-A3JHG_1_1107_7195_18281" 
#>  [567] "M00967_43_000000000-A3JHG_1_1107_7612_7844"  
#>  [568] "M00967_43_000000000-A3JHG_1_1107_7945_19118" 
#>  [569] "M00967_43_000000000-A3JHG_1_1107_7966_6427"  
#>  [570] "M00967_43_000000000-A3JHG_1_1107_8149_24158" 
#>  [571] "M00967_43_000000000-A3JHG_1_1107_8474_16938" 
#>  [572] "M00967_43_000000000-A3JHG_1_1107_9130_13455" 
#>  [573] "M00967_43_000000000-A3JHG_1_1107_9318_10158" 
#>  [574] "M00967_43_000000000-A3JHG_1_1108_10047_7859" 
#>  [575] "M00967_43_000000000-A3JHG_1_1108_10061_3691" 
#>  [576] "M00967_43_000000000-A3JHG_1_1108_10305_19494"
#>  [577] "M00967_43_000000000-A3JHG_1_1108_10511_13856"
#>  [578] "M00967_43_000000000-A3JHG_1_1108_10822_18112"
#>  [579] "M00967_43_000000000-A3JHG_1_1108_10911_24264"
#>  [580] "M00967_43_000000000-A3JHG_1_1108_11048_4561" 
#>  [581] "M00967_43_000000000-A3JHG_1_1108_11375_17687"
#>  [582] "M00967_43_000000000-A3JHG_1_1108_11391_5528" 
#>  [583] "M00967_43_000000000-A3JHG_1_1108_11559_26077"
#>  [584] "M00967_43_000000000-A3JHG_1_1108_11646_25892"
#>  [585] "M00967_43_000000000-A3JHG_1_1108_11757_26753"
#>  [586] "M00967_43_000000000-A3JHG_1_1108_12037_3053" 
#>  [587] "M00967_43_000000000-A3JHG_1_1108_12332_11456"
#>  [588] "M00967_43_000000000-A3JHG_1_1108_12507_2610" 
#>  [589] "M00967_43_000000000-A3JHG_1_1108_13402_18823"
#>  [590] "M00967_43_000000000-A3JHG_1_1108_13873_2429" 
#>  [591] "M00967_43_000000000-A3JHG_1_1108_14021_11937"
#>  [592] "M00967_43_000000000-A3JHG_1_1108_14282_18817"
#>  [593] "M00967_43_000000000-A3JHG_1_1108_14299_17220"
#>  [594] "M00967_43_000000000-A3JHG_1_1108_14370_23848"
#>  [595] "M00967_43_000000000-A3JHG_1_1108_14383_15110"
#>  [596] "M00967_43_000000000-A3JHG_1_1108_14399_23420"
#>  [597] "M00967_43_000000000-A3JHG_1_1108_14909_4503" 
#>  [598] "M00967_43_000000000-A3JHG_1_1108_15027_8794" 
#>  [599] "M00967_43_000000000-A3JHG_1_1108_15487_15318"
#>  [600] "M00967_43_000000000-A3JHG_1_1108_15711_24740"
#>  [601] "M00967_43_000000000-A3JHG_1_1108_15765_14436"
#>  [602] "M00967_43_000000000-A3JHG_1_1108_15844_3660" 
#>  [603] "M00967_43_000000000-A3JHG_1_1108_15903_3841" 
#>  [604] "M00967_43_000000000-A3JHG_1_1108_16241_15417"
#>  [605] "M00967_43_000000000-A3JHG_1_1108_16699_17018"
#>  [606] "M00967_43_000000000-A3JHG_1_1108_18362_23908"
#>  [607] "M00967_43_000000000-A3JHG_1_1108_18422_23262"
#>  [608] "M00967_43_000000000-A3JHG_1_1108_18644_12052"
#>  [609] "M00967_43_000000000-A3JHG_1_1108_18681_22857"
#>  [610] "M00967_43_000000000-A3JHG_1_1108_18791_11806"
#>  [611] "M00967_43_000000000-A3JHG_1_1108_18821_5778" 
#>  [612] "M00967_43_000000000-A3JHG_1_1108_18996_9436" 
#>  [613] "M00967_43_000000000-A3JHG_1_1108_19065_7362" 
#>  [614] "M00967_43_000000000-A3JHG_1_1108_19250_12721"
#>  [615] "M00967_43_000000000-A3JHG_1_1108_19698_21996"
#>  [616] "M00967_43_000000000-A3JHG_1_1108_19892_23978"
#>  [617] "M00967_43_000000000-A3JHG_1_1108_19974_13521"
#>  [618] "M00967_43_000000000-A3JHG_1_1108_20273_22925"
#>  [619] "M00967_43_000000000-A3JHG_1_1108_20614_22010"
#>  [620] "M00967_43_000000000-A3JHG_1_1108_20841_7623" 
#>  [621] "M00967_43_000000000-A3JHG_1_1108_21244_15970"
#>  [622] "M00967_43_000000000-A3JHG_1_1108_21553_18669"
#>  [623] "M00967_43_000000000-A3JHG_1_1108_21905_17218"
#>  [624] "M00967_43_000000000-A3JHG_1_1108_22541_9748" 
#>  [625] "M00967_43_000000000-A3JHG_1_1108_22866_20909"
#>  [626] "M00967_43_000000000-A3JHG_1_1108_23521_17609"
#>  [627] "M00967_43_000000000-A3JHG_1_1108_23570_6782" 
#>  [628] "M00967_43_000000000-A3JHG_1_1108_23963_16017"
#>  [629] "M00967_43_000000000-A3JHG_1_1108_23966_22243"
#>  [630] "M00967_43_000000000-A3JHG_1_1108_24102_4725" 
#>  [631] "M00967_43_000000000-A3JHG_1_1108_2439_18628" 
#>  [632] "M00967_43_000000000-A3JHG_1_1108_24697_25438"
#>  [633] "M00967_43_000000000-A3JHG_1_1108_25457_17518"
#>  [634] "M00967_43_000000000-A3JHG_1_1108_25479_20227"
#>  [635] "M00967_43_000000000-A3JHG_1_1108_26047_23516"
#>  [636] "M00967_43_000000000-A3JHG_1_1108_26418_17152"
#>  [637] "M00967_43_000000000-A3JHG_1_1108_26681_13284"
#>  [638] "M00967_43_000000000-A3JHG_1_1108_2714_17713" 
#>  [639] "M00967_43_000000000-A3JHG_1_1108_27163_22727"
#>  [640] "M00967_43_000000000-A3JHG_1_1108_27414_14298"
#>  [641] "M00967_43_000000000-A3JHG_1_1108_27464_20152"
#>  [642] "M00967_43_000000000-A3JHG_1_1108_27653_13933"
#>  [643] "M00967_43_000000000-A3JHG_1_1108_27927_21606"
#>  [644] "M00967_43_000000000-A3JHG_1_1108_28526_15940"
#>  [645] "M00967_43_000000000-A3JHG_1_1108_4844_10914" 
#>  [646] "M00967_43_000000000-A3JHG_1_1108_5894_23389" 
#>  [647] "M00967_43_000000000-A3JHG_1_1108_5951_24820" 
#>  [648] "M00967_43_000000000-A3JHG_1_1108_7641_20331" 
#>  [649] "M00967_43_000000000-A3JHG_1_1108_8714_18612" 
#>  [650] "M00967_43_000000000-A3JHG_1_1108_8770_12045" 
#>  [651] "M00967_43_000000000-A3JHG_1_1108_9199_25286" 
#>  [652] "M00967_43_000000000-A3JHG_1_1108_9239_22168" 
#>  [653] "M00967_43_000000000-A3JHG_1_1109_11007_22305"
#>  [654] "M00967_43_000000000-A3JHG_1_1109_11020_6991" 
#>  [655] "M00967_43_000000000-A3JHG_1_1109_11123_13961"
#>  [656] "M00967_43_000000000-A3JHG_1_1109_11295_11901"
#>  [657] "M00967_43_000000000-A3JHG_1_1109_11299_13989"
#>  [658] "M00967_43_000000000-A3JHG_1_1109_11326_21178"
#>  [659] "M00967_43_000000000-A3JHG_1_1109_12034_23177"
#>  [660] "M00967_43_000000000-A3JHG_1_1109_12129_4985" 
#>  [661] "M00967_43_000000000-A3JHG_1_1109_12151_14526"
#>  [662] "M00967_43_000000000-A3JHG_1_1109_12679_21321"
#>  [663] "M00967_43_000000000-A3JHG_1_1109_12715_19004"
#>  [664] "M00967_43_000000000-A3JHG_1_1109_12799_6921" 
#>  [665] "M00967_43_000000000-A3JHG_1_1109_13208_3584" 
#>  [666] "M00967_43_000000000-A3JHG_1_1109_13323_18613"
#>  [667] "M00967_43_000000000-A3JHG_1_1109_13330_21597"
#>  [668] "M00967_43_000000000-A3JHG_1_1109_13493_4507" 
#>  [669] "M00967_43_000000000-A3JHG_1_1109_13648_13972"
#>  [670] "M00967_43_000000000-A3JHG_1_1109_14359_12484"
#>  [671] "M00967_43_000000000-A3JHG_1_1109_14610_28830"
#>  [672] "M00967_43_000000000-A3JHG_1_1109_14735_26068"
#>  [673] "M00967_43_000000000-A3JHG_1_1109_15029_15883"
#>  [674] "M00967_43_000000000-A3JHG_1_1109_15598_4401" 
#>  [675] "M00967_43_000000000-A3JHG_1_1109_15821_9889" 
#>  [676] "M00967_43_000000000-A3JHG_1_1109_16375_7276" 
#>  [677] "M00967_43_000000000-A3JHG_1_1109_16714_23075"
#>  [678] "M00967_43_000000000-A3JHG_1_1109_17950_22790"
#>  [679] "M00967_43_000000000-A3JHG_1_1109_17953_9880" 
#>  [680] "M00967_43_000000000-A3JHG_1_1109_18126_7545" 
#>  [681] "M00967_43_000000000-A3JHG_1_1109_18164_26245"
#>  [682] "M00967_43_000000000-A3JHG_1_1109_18196_9880" 
#>  [683] "M00967_43_000000000-A3JHG_1_1109_18240_15644"
#>  [684] "M00967_43_000000000-A3JHG_1_1109_1885_14991" 
#>  [685] "M00967_43_000000000-A3JHG_1_1109_19298_22709"
#>  [686] "M00967_43_000000000-A3JHG_1_1109_19691_19130"
#>  [687] "M00967_43_000000000-A3JHG_1_1109_20401_24219"
#>  [688] "M00967_43_000000000-A3JHG_1_1109_20522_22061"
#>  [689] "M00967_43_000000000-A3JHG_1_1109_20716_17379"
#>  [690] "M00967_43_000000000-A3JHG_1_1109_21060_16840"
#>  [691] "M00967_43_000000000-A3JHG_1_1109_21458_15751"
#>  [692] "M00967_43_000000000-A3JHG_1_1109_21768_26571"
#>  [693] "M00967_43_000000000-A3JHG_1_1109_21900_16751"
#>  [694] "M00967_43_000000000-A3JHG_1_1109_22380_20532"
#>  [695] "M00967_43_000000000-A3JHG_1_1109_22483_18692"
#>  [696] "M00967_43_000000000-A3JHG_1_1109_22632_22311"
#>  [697] "M00967_43_000000000-A3JHG_1_1109_22795_14196"
#>  [698] "M00967_43_000000000-A3JHG_1_1109_23011_5924" 
#>  [699] "M00967_43_000000000-A3JHG_1_1109_23081_16781"
#>  [700] "M00967_43_000000000-A3JHG_1_1109_23085_25742"
#>  [701] "M00967_43_000000000-A3JHG_1_1109_23271_18418"
#>  [702] "M00967_43_000000000-A3JHG_1_1109_23907_9171" 
#>  [703] "M00967_43_000000000-A3JHG_1_1109_24068_9219" 
#>  [704] "M00967_43_000000000-A3JHG_1_1109_24109_8015" 
#>  [705] "M00967_43_000000000-A3JHG_1_1109_24246_18294"
#>  [706] "M00967_43_000000000-A3JHG_1_1109_24274_5733" 
#>  [707] "M00967_43_000000000-A3JHG_1_1109_24353_7102" 
#>  [708] "M00967_43_000000000-A3JHG_1_1109_24372_12546"
#>  [709] "M00967_43_000000000-A3JHG_1_1109_24493_23253"
#>  [710] "M00967_43_000000000-A3JHG_1_1109_24561_21388"
#>  [711] "M00967_43_000000000-A3JHG_1_1109_24712_5827" 
#>  [712] "M00967_43_000000000-A3JHG_1_1109_24850_14336"
#>  [713] "M00967_43_000000000-A3JHG_1_1109_26108_17578"
#>  [714] "M00967_43_000000000-A3JHG_1_1109_26230_18106"
#>  [715] "M00967_43_000000000-A3JHG_1_1109_26233_13220"
#>  [716] "M00967_43_000000000-A3JHG_1_1109_27552_22078"
#>  [717] "M00967_43_000000000-A3JHG_1_1109_2986_15638" 
#>  [718] "M00967_43_000000000-A3JHG_1_1109_3671_10080" 
#>  [719] "M00967_43_000000000-A3JHG_1_1109_4756_16948" 
#>  [720] "M00967_43_000000000-A3JHG_1_1109_5351_11203" 
#>  [721] "M00967_43_000000000-A3JHG_1_1109_5435_10528" 
#>  [722] "M00967_43_000000000-A3JHG_1_1109_6486_11603" 
#>  [723] "M00967_43_000000000-A3JHG_1_1109_6745_17470" 
#>  [724] "M00967_43_000000000-A3JHG_1_1109_6874_10395" 
#>  [725] "M00967_43_000000000-A3JHG_1_1109_7257_15397" 
#>  [726] "M00967_43_000000000-A3JHG_1_1109_7699_18608" 
#>  [727] "M00967_43_000000000-A3JHG_1_1109_7971_24259" 
#>  [728] "M00967_43_000000000-A3JHG_1_1109_7996_16802" 
#>  [729] "M00967_43_000000000-A3JHG_1_1109_8014_5010"  
#>  [730] "M00967_43_000000000-A3JHG_1_1109_8089_23582" 
#>  [731] "M00967_43_000000000-A3JHG_1_1109_8182_19921" 
#>  [732] "M00967_43_000000000-A3JHG_1_1109_8621_15527" 
#>  [733] "M00967_43_000000000-A3JHG_1_1109_8693_13422" 
#>  [734] "M00967_43_000000000-A3JHG_1_1109_8799_10155" 
#>  [735] "M00967_43_000000000-A3JHG_1_1109_8911_22008" 
#>  [736] "M00967_43_000000000-A3JHG_1_1109_8914_20463" 
#>  [737] "M00967_43_000000000-A3JHG_1_1109_8925_8850"  
#>  [738] "M00967_43_000000000-A3JHG_1_1109_9608_10994" 
#>  [739] "M00967_43_000000000-A3JHG_1_1109_9629_6769"  
#>  [740] "M00967_43_000000000-A3JHG_1_1109_9752_18363" 
#>  [741] "M00967_43_000000000-A3JHG_1_1110_10006_21831"
#>  [742] "M00967_43_000000000-A3JHG_1_1110_10033_25797"
#>  [743] "M00967_43_000000000-A3JHG_1_1110_10169_13508"
#>  [744] "M00967_43_000000000-A3JHG_1_1110_10693_14739"
#>  [745] "M00967_43_000000000-A3JHG_1_1110_10946_20858"
#>  [746] "M00967_43_000000000-A3JHG_1_1110_11254_14780"
#>  [747] "M00967_43_000000000-A3JHG_1_1110_11611_23800"
#>  [748] "M00967_43_000000000-A3JHG_1_1110_11891_7481" 
#>  [749] "M00967_43_000000000-A3JHG_1_1110_12199_12204"
#>  [750] "M00967_43_000000000-A3JHG_1_1110_12295_19699"
#>  [751] "M00967_43_000000000-A3JHG_1_1110_12737_15531"
#>  [752] "M00967_43_000000000-A3JHG_1_1110_12738_11413"
#>  [753] "M00967_43_000000000-A3JHG_1_1110_12935_6977" 
#>  [754] "M00967_43_000000000-A3JHG_1_1110_13033_21990"
#>  [755] "M00967_43_000000000-A3JHG_1_1110_13281_5144" 
#>  [756] "M00967_43_000000000-A3JHG_1_1110_13298_3746" 
#>  [757] "M00967_43_000000000-A3JHG_1_1110_13330_10638"
#>  [758] "M00967_43_000000000-A3JHG_1_1110_13719_21639"
#>  [759] "M00967_43_000000000-A3JHG_1_1110_13863_21427"
#>  [760] "M00967_43_000000000-A3JHG_1_1110_13879_1994" 
#>  [761] "M00967_43_000000000-A3JHG_1_1110_13931_12506"
#>  [762] "M00967_43_000000000-A3JHG_1_1110_13950_25143"
#>  [763] "M00967_43_000000000-A3JHG_1_1110_13999_26329"
#>  [764] "M00967_43_000000000-A3JHG_1_1110_14183_13016"
#>  [765] "M00967_43_000000000-A3JHG_1_1110_14278_11965"
#>  [766] "M00967_43_000000000-A3JHG_1_1110_14443_28339"
#>  [767] "M00967_43_000000000-A3JHG_1_1110_14465_21425"
#>  [768] "M00967_43_000000000-A3JHG_1_1110_14609_6975" 
#>  [769] "M00967_43_000000000-A3JHG_1_1110_14688_5568" 
#>  [770] "M00967_43_000000000-A3JHG_1_1110_14717_26485"
#>  [771] "M00967_43_000000000-A3JHG_1_1110_14814_10944"
#>  [772] "M00967_43_000000000-A3JHG_1_1110_14855_8462" 
#>  [773] "M00967_43_000000000-A3JHG_1_1110_15180_17571"
#>  [774] "M00967_43_000000000-A3JHG_1_1110_15277_6482" 
#>  [775] "M00967_43_000000000-A3JHG_1_1110_15335_17802"
#>  [776] "M00967_43_000000000-A3JHG_1_1110_15399_7509" 
#>  [777] "M00967_43_000000000-A3JHG_1_1110_15463_23056"
#>  [778] "M00967_43_000000000-A3JHG_1_1110_15641_10799"
#>  [779] "M00967_43_000000000-A3JHG_1_1110_15744_23114"
#>  [780] "M00967_43_000000000-A3JHG_1_1110_15833_26700"
#>  [781] "M00967_43_000000000-A3JHG_1_1110_15855_19234"
#>  [782] "M00967_43_000000000-A3JHG_1_1110_16116_28479"
#>  [783] "M00967_43_000000000-A3JHG_1_1110_16292_27414"
#>  [784] "M00967_43_000000000-A3JHG_1_1110_16691_26561"
#>  [785] "M00967_43_000000000-A3JHG_1_1110_16853_2840" 
#>  [786] "M00967_43_000000000-A3JHG_1_1110_17014_10530"
#>  [787] "M00967_43_000000000-A3JHG_1_1110_17088_10563"
#>  [788] "M00967_43_000000000-A3JHG_1_1110_17124_3344" 
#>  [789] "M00967_43_000000000-A3JHG_1_1110_17157_4300" 
#>  [790] "M00967_43_000000000-A3JHG_1_1110_17220_22884"
#>  [791] "M00967_43_000000000-A3JHG_1_1110_17368_18443"
#>  [792] "M00967_43_000000000-A3JHG_1_1110_17627_6261" 
#>  [793] "M00967_43_000000000-A3JHG_1_1110_17653_8119" 
#>  [794] "M00967_43_000000000-A3JHG_1_1110_17713_26784"
#>  [795] "M00967_43_000000000-A3JHG_1_1110_17805_7443" 
#>  [796] "M00967_43_000000000-A3JHG_1_1110_17818_21151"
#>  [797] "M00967_43_000000000-A3JHG_1_1110_17890_15072"
#>  [798] "M00967_43_000000000-A3JHG_1_1110_18108_15259"
#>  [799] "M00967_43_000000000-A3JHG_1_1110_18205_23899"
#>  [800] "M00967_43_000000000-A3JHG_1_1110_18313_3175" 
#>  [801] "M00967_43_000000000-A3JHG_1_1110_18344_10035"
#>  [802] "M00967_43_000000000-A3JHG_1_1110_18711_9702" 
#>  [803] "M00967_43_000000000-A3JHG_1_1110_18810_18264"
#>  [804] "M00967_43_000000000-A3JHG_1_1110_18818_19657"
#>  [805] "M00967_43_000000000-A3JHG_1_1110_18844_27717"
#>  [806] "M00967_43_000000000-A3JHG_1_1110_18950_6153" 
#>  [807] "M00967_43_000000000-A3JHG_1_1110_19425_20093"
#>  [808] "M00967_43_000000000-A3JHG_1_1110_19644_17655"
#>  [809] "M00967_43_000000000-A3JHG_1_1110_19955_4235" 
#>  [810] "M00967_43_000000000-A3JHG_1_1110_20049_16760"
#>  [811] "M00967_43_000000000-A3JHG_1_1110_20083_10445"
#>  [812] "M00967_43_000000000-A3JHG_1_1110_20095_8247" 
#>  [813] "M00967_43_000000000-A3JHG_1_1110_20540_12162"
#>  [814] "M00967_43_000000000-A3JHG_1_1110_20810_3468" 
#>  [815] "M00967_43_000000000-A3JHG_1_1110_20816_14150"
#>  [816] "M00967_43_000000000-A3JHG_1_1110_21383_22868"
#>  [817] "M00967_43_000000000-A3JHG_1_1110_21538_17471"
#>  [818] "M00967_43_000000000-A3JHG_1_1110_21542_3012" 
#>  [819] "M00967_43_000000000-A3JHG_1_1110_21619_15325"
#>  [820] "M00967_43_000000000-A3JHG_1_1110_21712_25652"
#>  [821] "M00967_43_000000000-A3JHG_1_1110_21794_19898"
#>  [822] "M00967_43_000000000-A3JHG_1_1110_21850_8958" 
#>  [823] "M00967_43_000000000-A3JHG_1_1110_22113_16232"
#>  [824] "M00967_43_000000000-A3JHG_1_1110_22175_11297"
#>  [825] "M00967_43_000000000-A3JHG_1_1110_22492_13115"
#>  [826] "M00967_43_000000000-A3JHG_1_1110_22982_25824"
#>  [827] "M00967_43_000000000-A3JHG_1_1110_23003_13541"
#>  [828] "M00967_43_000000000-A3JHG_1_1110_23557_12709"
#>  [829] "M00967_43_000000000-A3JHG_1_1110_24025_24881"
#>  [830] "M00967_43_000000000-A3JHG_1_1110_24287_19984"
#>  [831] "M00967_43_000000000-A3JHG_1_1110_24315_7382" 
#>  [832] "M00967_43_000000000-A3JHG_1_1110_24996_20556"
#>  [833] "M00967_43_000000000-A3JHG_1_1110_25043_9903" 
#>  [834] "M00967_43_000000000-A3JHG_1_1110_25186_13251"
#>  [835] "M00967_43_000000000-A3JHG_1_1110_25219_18794"
#>  [836] "M00967_43_000000000-A3JHG_1_1110_2564_19006" 
#>  [837] "M00967_43_000000000-A3JHG_1_1110_26069_24105"
#>  [838] "M00967_43_000000000-A3JHG_1_1110_26183_12645"
#>  [839] "M00967_43_000000000-A3JHG_1_1110_26301_16559"
#>  [840] "M00967_43_000000000-A3JHG_1_1110_26450_7049" 
#>  [841] "M00967_43_000000000-A3JHG_1_1110_26616_17589"
#>  [842] "M00967_43_000000000-A3JHG_1_1110_26933_11956"
#>  [843] "M00967_43_000000000-A3JHG_1_1110_26976_19838"
#>  [844] "M00967_43_000000000-A3JHG_1_1110_27186_10400"
#>  [845] "M00967_43_000000000-A3JHG_1_1110_27249_9297" 
#>  [846] "M00967_43_000000000-A3JHG_1_1110_27262_22663"
#>  [847] "M00967_43_000000000-A3JHG_1_1110_27306_15292"
#>  [848] "M00967_43_000000000-A3JHG_1_1110_27797_16769"
#>  [849] "M00967_43_000000000-A3JHG_1_1110_28012_17813"
#>  [850] "M00967_43_000000000-A3JHG_1_1110_28129_18678"
#>  [851] "M00967_43_000000000-A3JHG_1_1110_29167_14391"
#>  [852] "M00967_43_000000000-A3JHG_1_1110_3587_9347"  
#>  [853] "M00967_43_000000000-A3JHG_1_1110_3713_15410" 
#>  [854] "M00967_43_000000000-A3JHG_1_1110_4126_16552" 
#>  [855] "M00967_43_000000000-A3JHG_1_1110_4294_12570" 
#>  [856] "M00967_43_000000000-A3JHG_1_1110_4449_10146" 
#>  [857] "M00967_43_000000000-A3JHG_1_1110_4549_14740" 
#>  [858] "M00967_43_000000000-A3JHG_1_1110_4801_12425" 
#>  [859] "M00967_43_000000000-A3JHG_1_1110_4925_16493" 
#>  [860] "M00967_43_000000000-A3JHG_1_1110_5315_13833" 
#>  [861] "M00967_43_000000000-A3JHG_1_1110_5736_9217"  
#>  [862] "M00967_43_000000000-A3JHG_1_1110_5826_22164" 
#>  [863] "M00967_43_000000000-A3JHG_1_1110_6233_13348" 
#>  [864] "M00967_43_000000000-A3JHG_1_1110_6766_19046" 
#>  [865] "M00967_43_000000000-A3JHG_1_1110_6952_6984"  
#>  [866] "M00967_43_000000000-A3JHG_1_1110_7275_13443" 
#>  [867] "M00967_43_000000000-A3JHG_1_1110_7426_18774" 
#>  [868] "M00967_43_000000000-A3JHG_1_1110_7693_5256"  
#>  [869] "M00967_43_000000000-A3JHG_1_1110_7742_10390" 
#>  [870] "M00967_43_000000000-A3JHG_1_1110_7766_4240"  
#>  [871] "M00967_43_000000000-A3JHG_1_1110_7982_9474"  
#>  [872] "M00967_43_000000000-A3JHG_1_1110_8324_8698"  
#>  [873] "M00967_43_000000000-A3JHG_1_1110_8577_13143" 
#>  [874] "M00967_43_000000000-A3JHG_1_1110_8805_12834" 
#>  [875] "M00967_43_000000000-A3JHG_1_1110_9431_18167" 
#>  [876] "M00967_43_000000000-A3JHG_1_1110_9902_7407"  
#>  [877] "M00967_43_000000000-A3JHG_1_1111_10018_24068"
#>  [878] "M00967_43_000000000-A3JHG_1_1111_10106_14274"
#>  [879] "M00967_43_000000000-A3JHG_1_1111_10387_25149"
#>  [880] "M00967_43_000000000-A3JHG_1_1111_10475_27173"
#>  [881] "M00967_43_000000000-A3JHG_1_1111_10939_17099"
#>  [882] "M00967_43_000000000-A3JHG_1_1111_11315_12221"
#>  [883] "M00967_43_000000000-A3JHG_1_1111_11386_7041" 
#>  [884] "M00967_43_000000000-A3JHG_1_1111_11403_5737" 
#>  [885] "M00967_43_000000000-A3JHG_1_1111_11535_25208"
#>  [886] "M00967_43_000000000-A3JHG_1_1111_11611_14287"
#>  [887] "M00967_43_000000000-A3JHG_1_1111_11612_9640" 
#>  [888] "M00967_43_000000000-A3JHG_1_1111_11925_5550" 
#>  [889] "M00967_43_000000000-A3JHG_1_1111_12028_24705"
#>  [890] "M00967_43_000000000-A3JHG_1_1111_12107_16431"
#>  [891] "M00967_43_000000000-A3JHG_1_1111_12249_8143" 
#>  [892] "M00967_43_000000000-A3JHG_1_1111_12315_7486" 
#>  [893] "M00967_43_000000000-A3JHG_1_1111_12881_9910" 
#>  [894] "M00967_43_000000000-A3JHG_1_1111_12899_5675" 
#>  [895] "M00967_43_000000000-A3JHG_1_1111_13108_10829"
#>  [896] "M00967_43_000000000-A3JHG_1_1111_13451_28512"
#>  [897] "M00967_43_000000000-A3JHG_1_1111_13576_20790"
#>  [898] "M00967_43_000000000-A3JHG_1_1111_13586_2804" 
#>  [899] "M00967_43_000000000-A3JHG_1_1111_13633_13115"
#>  [900] "M00967_43_000000000-A3JHG_1_1111_13772_8177" 
#>  [901] "M00967_43_000000000-A3JHG_1_1111_13795_8174" 
#>  [902] "M00967_43_000000000-A3JHG_1_1111_13965_8116" 
#>  [903] "M00967_43_000000000-A3JHG_1_1111_14398_12314"
#>  [904] "M00967_43_000000000-A3JHG_1_1111_14408_20028"
#>  [905] "M00967_43_000000000-A3JHG_1_1111_15129_15739"
#>  [906] "M00967_43_000000000-A3JHG_1_1111_15247_9017" 
#>  [907] "M00967_43_000000000-A3JHG_1_1111_15312_18978"
#>  [908] "M00967_43_000000000-A3JHG_1_1111_15521_13523"
#>  [909] "M00967_43_000000000-A3JHG_1_1111_15637_23152"
#>  [910] "M00967_43_000000000-A3JHG_1_1111_16173_17418"
#>  [911] "M00967_43_000000000-A3JHG_1_1111_16200_26616"
#>  [912] "M00967_43_000000000-A3JHG_1_1111_16329_7024" 
#>  [913] "M00967_43_000000000-A3JHG_1_1111_16971_16405"
#>  [914] "M00967_43_000000000-A3JHG_1_1111_17002_7999" 
#>  [915] "M00967_43_000000000-A3JHG_1_1111_17009_11320"
#>  [916] "M00967_43_000000000-A3JHG_1_1111_17201_27600"
#>  [917] "M00967_43_000000000-A3JHG_1_1111_17207_6739" 
#>  [918] "M00967_43_000000000-A3JHG_1_1111_17550_17269"
#>  [919] "M00967_43_000000000-A3JHG_1_1111_17819_12937"
#>  [920] "M00967_43_000000000-A3JHG_1_1111_17856_12278"
#>  [921] "M00967_43_000000000-A3JHG_1_1111_18224_14666"
#>  [922] "M00967_43_000000000-A3JHG_1_1111_18784_23504"
#>  [923] "M00967_43_000000000-A3JHG_1_1111_18825_20030"
#>  [924] "M00967_43_000000000-A3JHG_1_1111_18978_23806"
#>  [925] "M00967_43_000000000-A3JHG_1_1111_19186_22787"
#>  [926] "M00967_43_000000000-A3JHG_1_1111_19384_17999"
#>  [927] "M00967_43_000000000-A3JHG_1_1111_19389_20281"
#>  [928] "M00967_43_000000000-A3JHG_1_1111_19645_10272"
#>  [929] "M00967_43_000000000-A3JHG_1_1111_20322_27077"
#>  [930] "M00967_43_000000000-A3JHG_1_1111_20588_10279"
#>  [931] "M00967_43_000000000-A3JHG_1_1111_20933_6700" 
#>  [932] "M00967_43_000000000-A3JHG_1_1111_21505_5523" 
#>  [933] "M00967_43_000000000-A3JHG_1_1111_21761_7263" 
#>  [934] "M00967_43_000000000-A3JHG_1_1111_21774_11978"
#>  [935] "M00967_43_000000000-A3JHG_1_1111_21999_17538"
#>  [936] "M00967_43_000000000-A3JHG_1_1111_22362_24215"
#>  [937] "M00967_43_000000000-A3JHG_1_1111_22910_7741" 
#>  [938] "M00967_43_000000000-A3JHG_1_1111_23503_22763"
#>  [939] "M00967_43_000000000-A3JHG_1_1111_23754_16069"
#>  [940] "M00967_43_000000000-A3JHG_1_1111_23893_19091"
#>  [941] "M00967_43_000000000-A3JHG_1_1111_24289_9499" 
#>  [942] "M00967_43_000000000-A3JHG_1_1111_24333_20558"
#>  [943] "M00967_43_000000000-A3JHG_1_1111_24385_9596" 
#>  [944] "M00967_43_000000000-A3JHG_1_1111_24465_10715"
#>  [945] "M00967_43_000000000-A3JHG_1_1111_24650_7997" 
#>  [946] "M00967_43_000000000-A3JHG_1_1111_24706_13976"
#>  [947] "M00967_43_000000000-A3JHG_1_1111_25166_7929" 
#>  [948] "M00967_43_000000000-A3JHG_1_1111_25344_22139"
#>  [949] "M00967_43_000000000-A3JHG_1_1111_25442_13485"
#>  [950] "M00967_43_000000000-A3JHG_1_1111_25638_13115"
#>  [951] "M00967_43_000000000-A3JHG_1_1111_25796_21952"
#>  [952] "M00967_43_000000000-A3JHG_1_1111_26760_12620"
#>  [953] "M00967_43_000000000-A3JHG_1_1111_29007_18395"
#>  [954] "M00967_43_000000000-A3JHG_1_1111_3610_17990" 
#>  [955] "M00967_43_000000000-A3JHG_1_1111_4607_20709" 
#>  [956] "M00967_43_000000000-A3JHG_1_1111_4667_8023"  
#>  [957] "M00967_43_000000000-A3JHG_1_1111_4823_22009" 
#>  [958] "M00967_43_000000000-A3JHG_1_1111_4893_9943"  
#>  [959] "M00967_43_000000000-A3JHG_1_1111_4974_8377"  
#>  [960] "M00967_43_000000000-A3JHG_1_1111_5281_15326" 
#>  [961] "M00967_43_000000000-A3JHG_1_1111_5848_24084" 
#>  [962] "M00967_43_000000000-A3JHG_1_1111_5994_19726" 
#>  [963] "M00967_43_000000000-A3JHG_1_1111_6250_5987"  
#>  [964] "M00967_43_000000000-A3JHG_1_1111_6317_24602" 
#>  [965] "M00967_43_000000000-A3JHG_1_1111_6871_23431" 
#>  [966] "M00967_43_000000000-A3JHG_1_1111_7157_8555"  
#>  [967] "M00967_43_000000000-A3JHG_1_1111_7963_16361" 
#>  [968] "M00967_43_000000000-A3JHG_1_1111_8141_21233" 
#>  [969] "M00967_43_000000000-A3JHG_1_1111_8276_5334"  
#>  [970] "M00967_43_000000000-A3JHG_1_1111_8296_18650" 
#>  [971] "M00967_43_000000000-A3JHG_1_1111_8373_10700" 
#>  [972] "M00967_43_000000000-A3JHG_1_1111_8697_7063"  
#>  [973] "M00967_43_000000000-A3JHG_1_1111_8849_21176" 
#>  [974] "M00967_43_000000000-A3JHG_1_1111_9294_7319"  
#>  [975] "M00967_43_000000000-A3JHG_1_1111_9296_21078" 
#>  [976] "M00967_43_000000000-A3JHG_1_1111_9420_7291"  
#>  [977] "M00967_43_000000000-A3JHG_1_1112_10024_21424"
#>  [978] "M00967_43_000000000-A3JHG_1_1112_10073_3247" 
#>  [979] "M00967_43_000000000-A3JHG_1_1112_10124_4953" 
#>  [980] "M00967_43_000000000-A3JHG_1_1112_10190_23278"
#>  [981] "M00967_43_000000000-A3JHG_1_1112_10194_20788"
#>  [982] "M00967_43_000000000-A3JHG_1_1112_10207_16135"
#>  [983] "M00967_43_000000000-A3JHG_1_1112_10447_18633"
#>  [984] "M00967_43_000000000-A3JHG_1_1112_10568_25353"
#>  [985] "M00967_43_000000000-A3JHG_1_1112_10641_14749"
#>  [986] "M00967_43_000000000-A3JHG_1_1112_10723_2426" 
#>  [987] "M00967_43_000000000-A3JHG_1_1112_11092_24430"
#>  [988] "M00967_43_000000000-A3JHG_1_1112_11192_11762"
#>  [989] "M00967_43_000000000-A3JHG_1_1112_11343_21537"
#>  [990] "M00967_43_000000000-A3JHG_1_1112_11370_15670"
#>  [991] "M00967_43_000000000-A3JHG_1_1112_11506_10347"
#>  [992] "M00967_43_000000000-A3JHG_1_1112_11788_5032" 
#>  [993] "M00967_43_000000000-A3JHG_1_1112_12079_14589"
#>  [994] "M00967_43_000000000-A3JHG_1_1112_12163_5660" 
#>  [995] "M00967_43_000000000-A3JHG_1_1112_12715_8487" 
#>  [996] "M00967_43_000000000-A3JHG_1_1112_12814_18299"
#>  [997] "M00967_43_000000000-A3JHG_1_1112_12844_3769" 
#>  [998] "M00967_43_000000000-A3JHG_1_1112_12844_6082" 
#>  [999] "M00967_43_000000000-A3JHG_1_1112_12887_9342" 
#> [1000] "M00967_43_000000000-A3JHG_1_1112_13142_18436"
#> [1001] "M00967_43_000000000-A3JHG_1_1112_13372_27170"
#> [1002] "M00967_43_000000000-A3JHG_1_1112_13563_20676"
#> [1003] "M00967_43_000000000-A3JHG_1_1112_13632_26764"
#> [1004] "M00967_43_000000000-A3JHG_1_1112_13976_5900" 
#> [1005] "M00967_43_000000000-A3JHG_1_1112_14080_4232" 
#> [1006] "M00967_43_000000000-A3JHG_1_1112_14201_3200" 
#> [1007] "M00967_43_000000000-A3JHG_1_1112_14508_26097"
#> [1008] "M00967_43_000000000-A3JHG_1_1112_14780_8078" 
#> [1009] "M00967_43_000000000-A3JHG_1_1112_14865_2720" 
#> [1010] "M00967_43_000000000-A3JHG_1_1112_15017_19202"
#> [1011] "M00967_43_000000000-A3JHG_1_1112_15029_25626"
#> [1012] "M00967_43_000000000-A3JHG_1_1112_15323_28130"
#> [1013] "M00967_43_000000000-A3JHG_1_1112_15603_7816" 
#> [1014] "M00967_43_000000000-A3JHG_1_1112_16159_23387"
#> [1015] "M00967_43_000000000-A3JHG_1_1112_16230_18061"
#> [1016] "M00967_43_000000000-A3JHG_1_1112_16388_28524"
#> [1017] "M00967_43_000000000-A3JHG_1_1112_16577_22524"
#> [1018] "M00967_43_000000000-A3JHG_1_1112_16673_28539"
#> [1019] "M00967_43_000000000-A3JHG_1_1112_16923_21328"
#> [1020] "M00967_43_000000000-A3JHG_1_1112_17696_15281"
#> [1021] "M00967_43_000000000-A3JHG_1_1112_17779_24453"
#> [1022] "M00967_43_000000000-A3JHG_1_1112_17857_20330"
#> [1023] "M00967_43_000000000-A3JHG_1_1112_17916_7193" 
#> [1024] "M00967_43_000000000-A3JHG_1_1112_17937_22981"
#> [1025] "M00967_43_000000000-A3JHG_1_1112_18213_15560"
#> [1026] "M00967_43_000000000-A3JHG_1_1112_18294_27835"
#> [1027] "M00967_43_000000000-A3JHG_1_1112_18378_12213"
#> [1028] "M00967_43_000000000-A3JHG_1_1112_18411_17052"
#> [1029] "M00967_43_000000000-A3JHG_1_1112_18485_9669" 
#> [1030] "M00967_43_000000000-A3JHG_1_1112_18506_25911"
#> [1031] "M00967_43_000000000-A3JHG_1_1112_18848_19266"
#> [1032] "M00967_43_000000000-A3JHG_1_1112_18981_20025"
#> [1033] "M00967_43_000000000-A3JHG_1_1112_19287_12246"
#> [1034] "M00967_43_000000000-A3JHG_1_1112_19376_17765"
#> [1035] "M00967_43_000000000-A3JHG_1_1112_19666_21477"
#> [1036] "M00967_43_000000000-A3JHG_1_1112_20127_11122"
#> [1037] "M00967_43_000000000-A3JHG_1_1112_20297_3328" 
#> [1038] "M00967_43_000000000-A3JHG_1_1112_20414_23391"
#> [1039] "M00967_43_000000000-A3JHG_1_1112_20872_8285" 
#> [1040] "M00967_43_000000000-A3JHG_1_1112_20884_13058"
#> [1041] "M00967_43_000000000-A3JHG_1_1112_20902_11558"
#> [1042] "M00967_43_000000000-A3JHG_1_1112_20972_19668"
#> [1043] "M00967_43_000000000-A3JHG_1_1112_21256_12063"
#> [1044] "M00967_43_000000000-A3JHG_1_1112_21353_8415" 
#> [1045] "M00967_43_000000000-A3JHG_1_1112_21719_22520"
#> [1046] "M00967_43_000000000-A3JHG_1_1112_22172_7939" 
#> [1047] "M00967_43_000000000-A3JHG_1_1112_22274_20191"
#> [1048] "M00967_43_000000000-A3JHG_1_1112_22736_6885" 
#> [1049] "M00967_43_000000000-A3JHG_1_1112_22806_22503"
#> [1050] "M00967_43_000000000-A3JHG_1_1112_22961_17195"
#> [1051] "M00967_43_000000000-A3JHG_1_1112_23531_7580" 
#> [1052] "M00967_43_000000000-A3JHG_1_1112_23738_7093" 
#> [1053] "M00967_43_000000000-A3JHG_1_1112_24118_19062"
#> [1054] "M00967_43_000000000-A3JHG_1_1112_24233_22070"
#> [1055] "M00967_43_000000000-A3JHG_1_1112_24606_18511"
#> [1056] "M00967_43_000000000-A3JHG_1_1112_24714_23272"
#> [1057] "M00967_43_000000000-A3JHG_1_1112_24720_18347"
#> [1058] "M00967_43_000000000-A3JHG_1_1112_24926_20157"
#> [1059] "M00967_43_000000000-A3JHG_1_1112_24992_23873"
#> [1060] "M00967_43_000000000-A3JHG_1_1112_25018_15007"
#> [1061] "M00967_43_000000000-A3JHG_1_1112_25108_12490"
#> [1062] "M00967_43_000000000-A3JHG_1_1112_25301_17173"
#> [1063] "M00967_43_000000000-A3JHG_1_1112_25336_20130"
#> [1064] "M00967_43_000000000-A3JHG_1_1112_25571_18822"
#> [1065] "M00967_43_000000000-A3JHG_1_1112_25719_18946"
#> [1066] "M00967_43_000000000-A3JHG_1_1112_25875_17229"
#> [1067] "M00967_43_000000000-A3JHG_1_1112_25924_15532"
#> [1068] "M00967_43_000000000-A3JHG_1_1112_26033_17779"
#> [1069] "M00967_43_000000000-A3JHG_1_1112_27130_19681"
#> [1070] "M00967_43_000000000-A3JHG_1_1112_27414_10900"
#> [1071] "M00967_43_000000000-A3JHG_1_1112_27485_13416"
#> [1072] "M00967_43_000000000-A3JHG_1_1112_27925_17349"
#> [1073] "M00967_43_000000000-A3JHG_1_1112_27973_10692"
#> [1074] "M00967_43_000000000-A3JHG_1_1112_28046_12095"
#> [1075] "M00967_43_000000000-A3JHG_1_1112_28222_11771"
#> [1076] "M00967_43_000000000-A3JHG_1_1112_28806_17427"
#> [1077] "M00967_43_000000000-A3JHG_1_1112_29269_13288"
#> [1078] "M00967_43_000000000-A3JHG_1_1112_2989_20122" 
#> [1079] "M00967_43_000000000-A3JHG_1_1112_3066_12463" 
#> [1080] "M00967_43_000000000-A3JHG_1_1112_5118_11366" 
#> [1081] "M00967_43_000000000-A3JHG_1_1112_5319_20617" 
#> [1082] "M00967_43_000000000-A3JHG_1_1112_5497_24301" 
#> [1083] "M00967_43_000000000-A3JHG_1_1112_5526_15494" 
#> [1084] "M00967_43_000000000-A3JHG_1_1112_5981_8948"  
#> [1085] "M00967_43_000000000-A3JHG_1_1112_6577_24020" 
#> [1086] "M00967_43_000000000-A3JHG_1_1112_6751_17524" 
#> [1087] "M00967_43_000000000-A3JHG_1_1112_6862_18037" 
#> [1088] "M00967_43_000000000-A3JHG_1_1112_6951_12447" 
#> [1089] "M00967_43_000000000-A3JHG_1_1112_7075_7100"  
#> [1090] "M00967_43_000000000-A3JHG_1_1112_7560_26006" 
#> [1091] "M00967_43_000000000-A3JHG_1_1112_7981_6310"  
#> [1092] "M00967_43_000000000-A3JHG_1_1112_8674_21438" 
#> [1093] "M00967_43_000000000-A3JHG_1_1112_8821_13863" 
#> [1094] "M00967_43_000000000-A3JHG_1_1112_9141_21408" 
#> [1095] "M00967_43_000000000-A3JHG_1_1112_9301_13162" 
#> [1096] "M00967_43_000000000-A3JHG_1_1112_9306_24571" 
#> [1097] "M00967_43_000000000-A3JHG_1_1112_9419_4295"  
#> [1098] "M00967_43_000000000-A3JHG_1_1112_9549_5158"  
#> [1099] "M00967_43_000000000-A3JHG_1_1112_9624_20695" 
#> [1100] "M00967_43_000000000-A3JHG_1_1112_9636_6179"  
#> [1101] "M00967_43_000000000-A3JHG_1_1112_9943_15231" 
#> [1102] "M00967_43_000000000-A3JHG_1_1113_10248_22431"
#> [1103] "M00967_43_000000000-A3JHG_1_1113_10719_18379"
#> [1104] "M00967_43_000000000-A3JHG_1_1113_11232_14725"
#> [1105] "M00967_43_000000000-A3JHG_1_1113_11248_15992"
#> [1106] "M00967_43_000000000-A3JHG_1_1113_11281_2740" 
#> [1107] "M00967_43_000000000-A3JHG_1_1113_11294_24024"
#> [1108] "M00967_43_000000000-A3JHG_1_1113_11760_14772"
#> [1109] "M00967_43_000000000-A3JHG_1_1113_11975_6040" 
#> [1110] "M00967_43_000000000-A3JHG_1_1113_12104_26590"
#> [1111] "M00967_43_000000000-A3JHG_1_1113_12514_11568"
#> [1112] "M00967_43_000000000-A3JHG_1_1113_12649_18573"
#> [1113] "M00967_43_000000000-A3JHG_1_1113_12711_3318" 
#> [1114] "M00967_43_000000000-A3JHG_1_1113_12882_9143" 
#> [1115] "M00967_43_000000000-A3JHG_1_1113_13105_21336"
#> [1116] "M00967_43_000000000-A3JHG_1_1113_13342_26294"
#> [1117] "M00967_43_000000000-A3JHG_1_1113_13413_8142" 
#> [1118] "M00967_43_000000000-A3JHG_1_1113_13627_24399"
#> [1119] "M00967_43_000000000-A3JHG_1_1113_13777_28036"
#> [1120] "M00967_43_000000000-A3JHG_1_1113_13883_23943"
#> [1121] "M00967_43_000000000-A3JHG_1_1113_13971_21290"
#> [1122] "M00967_43_000000000-A3JHG_1_1113_13987_16790"
#> [1123] "M00967_43_000000000-A3JHG_1_1113_14056_6035" 
#> [1124] "M00967_43_000000000-A3JHG_1_1113_14224_15638"
#> [1125] "M00967_43_000000000-A3JHG_1_1113_14260_20208"
#> [1126] "M00967_43_000000000-A3JHG_1_1113_14338_24655"
#> [1127] "M00967_43_000000000-A3JHG_1_1113_14485_10905"
#> [1128] "M00967_43_000000000-A3JHG_1_1113_14620_25479"
#> [1129] "M00967_43_000000000-A3JHG_1_1113_15332_23431"
#> [1130] "M00967_43_000000000-A3JHG_1_1113_15543_24961"
#> [1131] "M00967_43_000000000-A3JHG_1_1113_15548_11732"
#> [1132] "M00967_43_000000000-A3JHG_1_1113_15701_3674" 
#> [1133] "M00967_43_000000000-A3JHG_1_1113_15937_19783"
#> [1134] "M00967_43_000000000-A3JHG_1_1113_15967_2074" 
#> [1135] "M00967_43_000000000-A3JHG_1_1113_16063_11536"
#> [1136] "M00967_43_000000000-A3JHG_1_1113_16174_11502"
#> [1137] "M00967_43_000000000-A3JHG_1_1113_16472_4147" 
#> [1138] "M00967_43_000000000-A3JHG_1_1113_16596_16424"
#> [1139] "M00967_43_000000000-A3JHG_1_1113_16913_23641"
#> [1140] "M00967_43_000000000-A3JHG_1_1113_17095_9759" 
#> [1141] "M00967_43_000000000-A3JHG_1_1113_17274_13189"
#> [1142] "M00967_43_000000000-A3JHG_1_1113_17526_24276"
#> [1143] "M00967_43_000000000-A3JHG_1_1113_17833_23641"
#> [1144] "M00967_43_000000000-A3JHG_1_1113_18037_24127"
#> [1145] "M00967_43_000000000-A3JHG_1_1113_18114_26077"
#> [1146] "M00967_43_000000000-A3JHG_1_1113_18370_24961"
#> [1147] "M00967_43_000000000-A3JHG_1_1113_18462_22624"
#> [1148] "M00967_43_000000000-A3JHG_1_1113_18521_9535" 
#> [1149] "M00967_43_000000000-A3JHG_1_1113_18548_4292" 
#> [1150] "M00967_43_000000000-A3JHG_1_1113_18644_24962"
#> [1151] "M00967_43_000000000-A3JHG_1_1113_18958_10702"
#> [1152] "M00967_43_000000000-A3JHG_1_1113_19158_8422" 
#> [1153] "M00967_43_000000000-A3JHG_1_1113_19332_23681"
#> [1154] "M00967_43_000000000-A3JHG_1_1113_19457_3875" 
#> [1155] "M00967_43_000000000-A3JHG_1_1113_19532_4899" 
#> [1156] "M00967_43_000000000-A3JHG_1_1113_19563_8007" 
#> [1157] "M00967_43_000000000-A3JHG_1_1113_19689_14531"
#> [1158] "M00967_43_000000000-A3JHG_1_1113_19795_24825"
#> [1159] "M00967_43_000000000-A3JHG_1_1113_19958_8814" 
#> [1160] "M00967_43_000000000-A3JHG_1_1113_20185_10282"
#> [1161] "M00967_43_000000000-A3JHG_1_1113_20274_18664"
#> [1162] "M00967_43_000000000-A3JHG_1_1113_20770_14973"
#> [1163] "M00967_43_000000000-A3JHG_1_1113_20973_15490"
#> [1164] "M00967_43_000000000-A3JHG_1_1113_2128_16290" 
#> [1165] "M00967_43_000000000-A3JHG_1_1113_22254_26020"
#> [1166] "M00967_43_000000000-A3JHG_1_1113_22462_8142" 
#> [1167] "M00967_43_000000000-A3JHG_1_1113_22863_10725"
#> [1168] "M00967_43_000000000-A3JHG_1_1113_23113_7414" 
#> [1169] "M00967_43_000000000-A3JHG_1_1113_23401_8348" 
#> [1170] "M00967_43_000000000-A3JHG_1_1113_23502_11596"
#> [1171] "M00967_43_000000000-A3JHG_1_1113_23768_16589"
#> [1172] "M00967_43_000000000-A3JHG_1_1113_23986_14034"
#> [1173] "M00967_43_000000000-A3JHG_1_1113_24017_25976"
#> [1174] "M00967_43_000000000-A3JHG_1_1113_24111_15349"
#> [1175] "M00967_43_000000000-A3JHG_1_1113_24501_12692"
#> [1176] "M00967_43_000000000-A3JHG_1_1113_25205_6894" 
#> [1177] "M00967_43_000000000-A3JHG_1_1113_25390_24192"
#> [1178] "M00967_43_000000000-A3JHG_1_1113_25479_14916"
#> [1179] "M00967_43_000000000-A3JHG_1_1113_25500_12155"
#> [1180] "M00967_43_000000000-A3JHG_1_1113_25843_7710" 
#> [1181] "M00967_43_000000000-A3JHG_1_1113_25861_7699" 
#> [1182] "M00967_43_000000000-A3JHG_1_1113_26249_6887" 
#> [1183] "M00967_43_000000000-A3JHG_1_1113_26510_16371"
#> [1184] "M00967_43_000000000-A3JHG_1_1113_26517_18281"
#> [1185] "M00967_43_000000000-A3JHG_1_1113_26736_16107"
#> [1186] "M00967_43_000000000-A3JHG_1_1113_26886_12231"
#> [1187] "M00967_43_000000000-A3JHG_1_1113_27215_10406"
#> [1188] "M00967_43_000000000-A3JHG_1_1113_27297_13091"
#> [1189] "M00967_43_000000000-A3JHG_1_1113_28094_17280"
#> [1190] "M00967_43_000000000-A3JHG_1_1113_28261_14681"
#> [1191] "M00967_43_000000000-A3JHG_1_1113_28574_18008"
#> [1192] "M00967_43_000000000-A3JHG_1_1113_28928_16770"
#> [1193] "M00967_43_000000000-A3JHG_1_1113_29188_13629"
#> [1194] "M00967_43_000000000-A3JHG_1_1113_4318_16941" 
#> [1195] "M00967_43_000000000-A3JHG_1_1113_4659_12375" 
#> [1196] "M00967_43_000000000-A3JHG_1_1113_5336_24219" 
#> [1197] "M00967_43_000000000-A3JHG_1_1113_5425_15020" 
#> [1198] "M00967_43_000000000-A3JHG_1_1113_6268_10349" 
#> [1199] "M00967_43_000000000-A3JHG_1_1113_6303_20947" 
#> [1200] "M00967_43_000000000-A3JHG_1_1113_6524_23905" 
#> [1201] "M00967_43_000000000-A3JHG_1_1113_6556_22579" 
#> [1202] "M00967_43_000000000-A3JHG_1_1113_6843_21123" 
#> [1203] "M00967_43_000000000-A3JHG_1_1113_7278_14313" 
#> [1204] "M00967_43_000000000-A3JHG_1_1113_7625_23288" 
#> [1205] "M00967_43_000000000-A3JHG_1_1113_8139_17897" 
#> [1206] "M00967_43_000000000-A3JHG_1_1113_8282_21839" 
#> [1207] "M00967_43_000000000-A3JHG_1_1113_8308_7055"  
#> [1208] "M00967_43_000000000-A3JHG_1_1113_8392_20172" 
#> [1209] "M00967_43_000000000-A3JHG_1_1113_9185_6542"  
#> [1210] "M00967_43_000000000-A3JHG_1_1113_9276_24879" 
#> [1211] "M00967_43_000000000-A3JHG_1_1113_9372_3141"  
#> [1212] "M00967_43_000000000-A3JHG_1_1113_9498_21548" 
#> [1213] "M00967_43_000000000-A3JHG_1_1114_10011_16636"
#> [1214] "M00967_43_000000000-A3JHG_1_1114_10055_21108"
#> [1215] "M00967_43_000000000-A3JHG_1_1114_10193_19783"
#> [1216] "M00967_43_000000000-A3JHG_1_1114_10373_6011" 
#> [1217] "M00967_43_000000000-A3JHG_1_1114_10520_9402" 
#> [1218] "M00967_43_000000000-A3JHG_1_1114_10573_19345"
#> [1219] "M00967_43_000000000-A3JHG_1_1114_10874_10296"
#> [1220] "M00967_43_000000000-A3JHG_1_1114_11189_15616"
#> [1221] "M00967_43_000000000-A3JHG_1_1114_11237_4361" 
#> [1222] "M00967_43_000000000-A3JHG_1_1114_11769_3401" 
#> [1223] "M00967_43_000000000-A3JHG_1_1114_11845_9406" 
#> [1224] "M00967_43_000000000-A3JHG_1_1114_12174_8389" 
#> [1225] "M00967_43_000000000-A3JHG_1_1114_12357_13871"
#> [1226] "M00967_43_000000000-A3JHG_1_1114_12471_18325"
#> [1227] "M00967_43_000000000-A3JHG_1_1114_12578_20022"
#> [1228] "M00967_43_000000000-A3JHG_1_1114_12644_14109"
#> [1229] "M00967_43_000000000-A3JHG_1_1114_12898_16926"
#> [1230] "M00967_43_000000000-A3JHG_1_1114_12965_2018" 
#> [1231] "M00967_43_000000000-A3JHG_1_1114_13158_9182" 
#> [1232] "M00967_43_000000000-A3JHG_1_1114_13201_9116" 
#> [1233] "M00967_43_000000000-A3JHG_1_1114_13258_26764"
#> [1234] "M00967_43_000000000-A3JHG_1_1114_13266_17804"
#> [1235] "M00967_43_000000000-A3JHG_1_1114_13389_13295"
#> [1236] "M00967_43_000000000-A3JHG_1_1114_13397_2211" 
#> [1237] "M00967_43_000000000-A3JHG_1_1114_13512_6898" 
#> [1238] "M00967_43_000000000-A3JHG_1_1114_13556_18457"
#> [1239] "M00967_43_000000000-A3JHG_1_1114_13941_17165"
#> [1240] "M00967_43_000000000-A3JHG_1_1114_13983_6574" 
#> [1241] "M00967_43_000000000-A3JHG_1_1114_14259_12215"
#> [1242] "M00967_43_000000000-A3JHG_1_1114_14280_26421"
#> [1243] "M00967_43_000000000-A3JHG_1_1114_14308_17024"
#> [1244] "M00967_43_000000000-A3JHG_1_1114_14431_2336" 
#> [1245] "M00967_43_000000000-A3JHG_1_1114_14691_17449"
#> [1246] "M00967_43_000000000-A3JHG_1_1114_14693_26365"
#> [1247] "M00967_43_000000000-A3JHG_1_1114_14716_28079"
#> [1248] "M00967_43_000000000-A3JHG_1_1114_14939_20476"
#> [1249] "M00967_43_000000000-A3JHG_1_1114_14958_19894"
#> [1250] "M00967_43_000000000-A3JHG_1_1114_14996_3562" 
#> [1251] "M00967_43_000000000-A3JHG_1_1114_15218_23608"
#> [1252] "M00967_43_000000000-A3JHG_1_1114_15285_11132"
#> [1253] "M00967_43_000000000-A3JHG_1_1114_15409_10721"
#> [1254] "M00967_43_000000000-A3JHG_1_1114_15498_19369"
#> [1255] "M00967_43_000000000-A3JHG_1_1114_15727_25995"
#> [1256] "M00967_43_000000000-A3JHG_1_1114_16061_15885"
#> [1257] "M00967_43_000000000-A3JHG_1_1114_16110_14718"
#> [1258] "M00967_43_000000000-A3JHG_1_1114_16584_12333"
#> [1259] "M00967_43_000000000-A3JHG_1_1114_16655_8470" 
#> [1260] "M00967_43_000000000-A3JHG_1_1114_16967_16702"
#> [1261] "M00967_43_000000000-A3JHG_1_1114_17730_26068"
#> [1262] "M00967_43_000000000-A3JHG_1_1114_17780_27643"
#> [1263] "M00967_43_000000000-A3JHG_1_1114_17898_13576"
#> [1264] "M00967_43_000000000-A3JHG_1_1114_18112_24383"
#> [1265] "M00967_43_000000000-A3JHG_1_1114_18273_11520"
#> [1266] "M00967_43_000000000-A3JHG_1_1114_18296_7709" 
#> [1267] "M00967_43_000000000-A3JHG_1_1114_18510_8776" 
#> [1268] "M00967_43_000000000-A3JHG_1_1114_1863_15569" 
#> [1269] "M00967_43_000000000-A3JHG_1_1114_18788_21055"
#> [1270] "M00967_43_000000000-A3JHG_1_1114_18793_4661" 
#> [1271] "M00967_43_000000000-A3JHG_1_1114_19154_8327" 
#> [1272] "M00967_43_000000000-A3JHG_1_1114_19506_16402"
#> [1273] "M00967_43_000000000-A3JHG_1_1114_19514_5000" 
#> [1274] "M00967_43_000000000-A3JHG_1_1114_19529_17837"
#> [1275] "M00967_43_000000000-A3JHG_1_1114_19698_19740"
#> [1276] "M00967_43_000000000-A3JHG_1_1114_19795_28223"
#> [1277] "M00967_43_000000000-A3JHG_1_1114_19875_9916" 
#> [1278] "M00967_43_000000000-A3JHG_1_1114_20218_10272"
#> [1279] "M00967_43_000000000-A3JHG_1_1114_20220_13470"
#> [1280] "M00967_43_000000000-A3JHG_1_1114_20440_13060"
#> [1281] "M00967_43_000000000-A3JHG_1_1114_21222_13060"
#> [1282] "M00967_43_000000000-A3JHG_1_1114_21515_19829"
#> [1283] "M00967_43_000000000-A3JHG_1_1114_21549_4010" 
#> [1284] "M00967_43_000000000-A3JHG_1_1114_21607_13167"
#> [1285] "M00967_43_000000000-A3JHG_1_1114_21630_21181"
#> [1286] "M00967_43_000000000-A3JHG_1_1114_21716_25171"
#> [1287] "M00967_43_000000000-A3JHG_1_1114_21924_11050"
#> [1288] "M00967_43_000000000-A3JHG_1_1114_22144_24942"
#> [1289] "M00967_43_000000000-A3JHG_1_1114_22494_17866"
#> [1290] "M00967_43_000000000-A3JHG_1_1114_22681_21130"
#> [1291] "M00967_43_000000000-A3JHG_1_1114_22706_5018" 
#> [1292] "M00967_43_000000000-A3JHG_1_1114_22804_12324"
#> [1293] "M00967_43_000000000-A3JHG_1_1114_22942_23318"
#> [1294] "M00967_43_000000000-A3JHG_1_1114_23025_24395"
#> [1295] "M00967_43_000000000-A3JHG_1_1114_23055_10307"
#> [1296] "M00967_43_000000000-A3JHG_1_1114_23407_6585" 
#> [1297] "M00967_43_000000000-A3JHG_1_1114_23518_14384"
#> [1298] "M00967_43_000000000-A3JHG_1_1114_23713_17193"
#> [1299] "M00967_43_000000000-A3JHG_1_1114_23881_16039"
#> [1300] "M00967_43_000000000-A3JHG_1_1114_23944_10783"
#> [1301] "M00967_43_000000000-A3JHG_1_1114_24088_5500" 
#> [1302] "M00967_43_000000000-A3JHG_1_1114_24430_23737"
#> [1303] "M00967_43_000000000-A3JHG_1_1114_24795_10657"
#> [1304] "M00967_43_000000000-A3JHG_1_1114_25025_23939"
#> [1305] "M00967_43_000000000-A3JHG_1_1114_25131_8732" 
#> [1306] "M00967_43_000000000-A3JHG_1_1114_25590_23883"
#> [1307] "M00967_43_000000000-A3JHG_1_1114_25962_10907"
#> [1308] "M00967_43_000000000-A3JHG_1_1114_26027_6706" 
#> [1309] "M00967_43_000000000-A3JHG_1_1114_26942_22032"
#> [1310] "M00967_43_000000000-A3JHG_1_1114_26968_12650"
#> [1311] "M00967_43_000000000-A3JHG_1_1114_27071_8397" 
#> [1312] "M00967_43_000000000-A3JHG_1_1114_27304_18058"
#> [1313] "M00967_43_000000000-A3JHG_1_1114_27471_20890"
#> [1314] "M00967_43_000000000-A3JHG_1_1114_27837_14578"
#> [1315] "M00967_43_000000000-A3JHG_1_1114_27907_14387"
#> [1316] "M00967_43_000000000-A3JHG_1_1114_28031_16211"
#> [1317] "M00967_43_000000000-A3JHG_1_1114_3131_16050" 
#> [1318] "M00967_43_000000000-A3JHG_1_1114_3366_17894" 
#> [1319] "M00967_43_000000000-A3JHG_1_1114_3483_21724" 
#> [1320] "M00967_43_000000000-A3JHG_1_1114_4370_15679" 
#> [1321] "M00967_43_000000000-A3JHG_1_1114_5207_13479" 
#> [1322] "M00967_43_000000000-A3JHG_1_1114_5681_17379" 
#> [1323] "M00967_43_000000000-A3JHG_1_1114_5788_19276" 
#> [1324] "M00967_43_000000000-A3JHG_1_1114_5902_6841"  
#> [1325] "M00967_43_000000000-A3JHG_1_1114_6141_22497" 
#> [1326] "M00967_43_000000000-A3JHG_1_1114_6229_5326"  
#> [1327] "M00967_43_000000000-A3JHG_1_1114_6495_13602" 
#> [1328] "M00967_43_000000000-A3JHG_1_1114_6688_5293"  
#> [1329] "M00967_43_000000000-A3JHG_1_1114_7103_21358" 
#> [1330] "M00967_43_000000000-A3JHG_1_1114_7341_25496" 
#> [1331] "M00967_43_000000000-A3JHG_1_1114_7390_15529" 
#> [1332] "M00967_43_000000000-A3JHG_1_1114_7614_20887" 
#> [1333] "M00967_43_000000000-A3JHG_1_1114_7676_25736" 
#> [1334] "M00967_43_000000000-A3JHG_1_1114_7920_25776" 
#> [1335] "M00967_43_000000000-A3JHG_1_1114_8059_18290" 
#> [1336] "M00967_43_000000000-A3JHG_1_1114_8142_5567"  
#> [1337] "M00967_43_000000000-A3JHG_1_1114_8377_12333" 
#> [1338] "M00967_43_000000000-A3JHG_1_1114_8492_23708" 
#> [1339] "M00967_43_000000000-A3JHG_1_1114_8547_14272" 
#> [1340] "M00967_43_000000000-A3JHG_1_1114_8550_6381"  
#> [1341] "M00967_43_000000000-A3JHG_1_1114_8603_11968" 
#> [1342] "M00967_43_000000000-A3JHG_1_1114_8628_21912" 
#> [1343] "M00967_43_000000000-A3JHG_1_1114_8676_16145" 
#> [1344] "M00967_43_000000000-A3JHG_1_1114_8829_4071"  
#> [1345] "M00967_43_000000000-A3JHG_1_1114_8934_14248" 
#> [1346] "M00967_43_000000000-A3JHG_1_1114_9024_20163" 
#> [1347] "M00967_43_000000000-A3JHG_1_1114_9210_8620"  
#> [1348] "M00967_43_000000000-A3JHG_1_1114_9365_4621"  
#> [1349] "M00967_43_000000000-A3JHG_1_1114_9368_6991"  
#> [1350] "M00967_43_000000000-A3JHG_1_1114_9494_11127" 
#> [1351] "M00967_43_000000000-A3JHG_1_1114_9494_6921"  
#> [1352] "M00967_43_000000000-A3JHG_1_2101_10084_22191"
#> [1353] "M00967_43_000000000-A3JHG_1_2101_10104_19430"
#> [1354] "M00967_43_000000000-A3JHG_1_2101_10680_16017"
#> [1355] "M00967_43_000000000-A3JHG_1_2101_10845_21395"
#> [1356] "M00967_43_000000000-A3JHG_1_2101_11015_9186" 
#> [1357] "M00967_43_000000000-A3JHG_1_2101_12238_5807" 
#> [1358] "M00967_43_000000000-A3JHG_1_2101_12255_19099"
#> [1359] "M00967_43_000000000-A3JHG_1_2101_12402_11394"
#> [1360] "M00967_43_000000000-A3JHG_1_2101_12505_13702"
#> [1361] "M00967_43_000000000-A3JHG_1_2101_13105_15429"
#> [1362] "M00967_43_000000000-A3JHG_1_2101_13130_17466"
#> [1363] "M00967_43_000000000-A3JHG_1_2101_13320_3436" 
#> [1364] "M00967_43_000000000-A3JHG_1_2101_13474_15202"
#> [1365] "M00967_43_000000000-A3JHG_1_2101_13595_15201"
#> [1366] "M00967_43_000000000-A3JHG_1_2101_13838_9027" 
#> [1367] "M00967_43_000000000-A3JHG_1_2101_13911_2481" 
#> [1368] "M00967_43_000000000-A3JHG_1_2101_13911_4901" 
#> [1369] "M00967_43_000000000-A3JHG_1_2101_14159_9619" 
#> [1370] "M00967_43_000000000-A3JHG_1_2101_14628_22758"
#> [1371] "M00967_43_000000000-A3JHG_1_2101_15190_13450"
#> [1372] "M00967_43_000000000-A3JHG_1_2101_16005_11965"
#> [1373] "M00967_43_000000000-A3JHG_1_2101_16339_18404"
#> [1374] "M00967_43_000000000-A3JHG_1_2101_16354_15467"
#> [1375] "M00967_43_000000000-A3JHG_1_2101_16474_12783"
#> [1376] "M00967_43_000000000-A3JHG_1_2101_16487_8614" 
#> [1377] "M00967_43_000000000-A3JHG_1_2101_16841_5683" 
#> [1378] "M00967_43_000000000-A3JHG_1_2101_16955_14005"
#> [1379] "M00967_43_000000000-A3JHG_1_2101_17219_12550"
#> [1380] "M00967_43_000000000-A3JHG_1_2101_17292_21391"
#> [1381] "M00967_43_000000000-A3JHG_1_2101_17405_14654"
#> [1382] "M00967_43_000000000-A3JHG_1_2101_17437_26807"
#> [1383] "M00967_43_000000000-A3JHG_1_2101_17475_7154" 
#> [1384] "M00967_43_000000000-A3JHG_1_2101_17535_7182" 
#> [1385] "M00967_43_000000000-A3JHG_1_2101_17971_14439"
#> [1386] "M00967_43_000000000-A3JHG_1_2101_20087_25441"
#> [1387] "M00967_43_000000000-A3JHG_1_2101_20145_22042"
#> [1388] "M00967_43_000000000-A3JHG_1_2101_20164_2274" 
#> [1389] "M00967_43_000000000-A3JHG_1_2101_20223_25753"
#> [1390] "M00967_43_000000000-A3JHG_1_2101_21575_17061"
#> [1391] "M00967_43_000000000-A3JHG_1_2101_21882_9957" 
#> [1392] "M00967_43_000000000-A3JHG_1_2101_21947_18735"
#> [1393] "M00967_43_000000000-A3JHG_1_2101_22000_19843"
#> [1394] "M00967_43_000000000-A3JHG_1_2101_22171_7080" 
#> [1395] "M00967_43_000000000-A3JHG_1_2101_22400_13416"
#> [1396] "M00967_43_000000000-A3JHG_1_2101_22555_16416"
#> [1397] "M00967_43_000000000-A3JHG_1_2101_22573_22946"
#> [1398] "M00967_43_000000000-A3JHG_1_2101_22942_18205"
#> [1399] "M00967_43_000000000-A3JHG_1_2101_22980_4780" 
#> [1400] "M00967_43_000000000-A3JHG_1_2101_23114_22755"
#> [1401] "M00967_43_000000000-A3JHG_1_2101_23950_13306"
#> [1402] "M00967_43_000000000-A3JHG_1_2101_24471_17776"
#> [1403] "M00967_43_000000000-A3JHG_1_2101_24821_18536"
#> [1404] "M00967_43_000000000-A3JHG_1_2101_24905_9024" 
#> [1405] "M00967_43_000000000-A3JHG_1_2101_25094_24210"
#> [1406] "M00967_43_000000000-A3JHG_1_2101_25211_5709" 
#> [1407] "M00967_43_000000000-A3JHG_1_2101_25412_7972" 
#> [1408] "M00967_43_000000000-A3JHG_1_2101_25710_23156"
#> [1409] "M00967_43_000000000-A3JHG_1_2101_25853_16524"
#> [1410] "M00967_43_000000000-A3JHG_1_2101_26265_15434"
#> [1411] "M00967_43_000000000-A3JHG_1_2101_26453_14067"
#> [1412] "M00967_43_000000000-A3JHG_1_2101_26706_7767" 
#> [1413] "M00967_43_000000000-A3JHG_1_2101_26784_17549"
#> [1414] "M00967_43_000000000-A3JHG_1_2101_27541_14282"
#> [1415] "M00967_43_000000000-A3JHG_1_2101_27710_20218"
#> [1416] "M00967_43_000000000-A3JHG_1_2101_27920_18472"
#> [1417] "M00967_43_000000000-A3JHG_1_2101_28033_9425" 
#> [1418] "M00967_43_000000000-A3JHG_1_2101_28391_18939"
#> [1419] "M00967_43_000000000-A3JHG_1_2101_3680_12611" 
#> [1420] "M00967_43_000000000-A3JHG_1_2101_4117_14210" 
#> [1421] "M00967_43_000000000-A3JHG_1_2101_4159_12177" 
#> [1422] "M00967_43_000000000-A3JHG_1_2101_4572_15105" 
#> [1423] "M00967_43_000000000-A3JHG_1_2101_4780_15092" 
#> [1424] "M00967_43_000000000-A3JHG_1_2101_6223_8694"  
#> [1425] "M00967_43_000000000-A3JHG_1_2101_6377_14009" 
#> [1426] "M00967_43_000000000-A3JHG_1_2101_6909_12719" 
#> [1427] "M00967_43_000000000-A3JHG_1_2101_7036_14239" 
#> [1428] "M00967_43_000000000-A3JHG_1_2101_7254_13419" 
#> [1429] "M00967_43_000000000-A3JHG_1_2101_7382_19629" 
#> [1430] "M00967_43_000000000-A3JHG_1_2101_7947_10616" 
#> [1431] "M00967_43_000000000-A3JHG_1_2101_8442_23370" 
#> [1432] "M00967_43_000000000-A3JHG_1_2101_9084_14502" 
#> [1433] "M00967_43_000000000-A3JHG_1_2101_9500_21244" 
#> [1434] "M00967_43_000000000-A3JHG_1_2101_9650_25321" 
#> [1435] "M00967_43_000000000-A3JHG_1_2102_10158_13191"
#> [1436] "M00967_43_000000000-A3JHG_1_2102_10617_12843"
#> [1437] "M00967_43_000000000-A3JHG_1_2102_10627_9330" 
#> [1438] "M00967_43_000000000-A3JHG_1_2102_10683_16772"
#> [1439] "M00967_43_000000000-A3JHG_1_2102_11266_17437"
#> [1440] "M00967_43_000000000-A3JHG_1_2102_11553_9071" 
#> [1441] "M00967_43_000000000-A3JHG_1_2102_11614_20083"
#> [1442] "M00967_43_000000000-A3JHG_1_2102_12098_12991"
#> [1443] "M00967_43_000000000-A3JHG_1_2102_12263_16584"
#> [1444] "M00967_43_000000000-A3JHG_1_2102_12357_13497"
#> [1445] "M00967_43_000000000-A3JHG_1_2102_12696_27241"
#> [1446] "M00967_43_000000000-A3JHG_1_2102_12732_5062" 
#> [1447] "M00967_43_000000000-A3JHG_1_2102_12782_8845" 
#> [1448] "M00967_43_000000000-A3JHG_1_2102_12937_4366" 
#> [1449] "M00967_43_000000000-A3JHG_1_2102_13864_6831" 
#> [1450] "M00967_43_000000000-A3JHG_1_2102_14252_8710" 
#> [1451] "M00967_43_000000000-A3JHG_1_2102_14512_24313"
#> [1452] "M00967_43_000000000-A3JHG_1_2102_14650_8724" 
#> [1453] "M00967_43_000000000-A3JHG_1_2102_15356_3054" 
#> [1454] "M00967_43_000000000-A3JHG_1_2102_15542_26928"
#> [1455] "M00967_43_000000000-A3JHG_1_2102_15674_19713"
#> [1456] "M00967_43_000000000-A3JHG_1_2102_15692_18713"
#> [1457] "M00967_43_000000000-A3JHG_1_2102_15728_25058"
#> [1458] "M00967_43_000000000-A3JHG_1_2102_15778_12420"
#> [1459] "M00967_43_000000000-A3JHG_1_2102_16184_16848"
#> [1460] "M00967_43_000000000-A3JHG_1_2102_16491_2474" 
#> [1461] "M00967_43_000000000-A3JHG_1_2102_16703_5011" 
#> [1462] "M00967_43_000000000-A3JHG_1_2102_16879_22327"
#> [1463] "M00967_43_000000000-A3JHG_1_2102_16924_2551" 
#> [1464] "M00967_43_000000000-A3JHG_1_2102_17019_10818"
#> [1465] "M00967_43_000000000-A3JHG_1_2102_17055_16238"
#> [1466] "M00967_43_000000000-A3JHG_1_2102_17280_15307"
#> [1467] "M00967_43_000000000-A3JHG_1_2102_17332_24887"
#> [1468] "M00967_43_000000000-A3JHG_1_2102_17415_24873"
#> [1469] "M00967_43_000000000-A3JHG_1_2102_17714_13657"
#> [1470] "M00967_43_000000000-A3JHG_1_2102_17791_24913"
#> [1471] "M00967_43_000000000-A3JHG_1_2102_18322_3392" 
#> [1472] "M00967_43_000000000-A3JHG_1_2102_18868_8410" 
#> [1473] "M00967_43_000000000-A3JHG_1_2102_19287_23968"
#> [1474] "M00967_43_000000000-A3JHG_1_2102_19327_26344"
#> [1475] "M00967_43_000000000-A3JHG_1_2102_21011_11727"
#> [1476] "M00967_43_000000000-A3JHG_1_2102_21189_15385"
#> [1477] "M00967_43_000000000-A3JHG_1_2102_21816_16599"
#> [1478] "M00967_43_000000000-A3JHG_1_2102_22092_15614"
#> [1479] "M00967_43_000000000-A3JHG_1_2102_22167_26446"
#> [1480] "M00967_43_000000000-A3JHG_1_2102_22381_10329"
#> [1481] "M00967_43_000000000-A3JHG_1_2102_22437_24206"
#> [1482] "M00967_43_000000000-A3JHG_1_2102_22604_4240" 
#> [1483] "M00967_43_000000000-A3JHG_1_2102_23398_14571"
#> [1484] "M00967_43_000000000-A3JHG_1_2102_24666_20800"
#> [1485] "M00967_43_000000000-A3JHG_1_2102_24708_11365"
#> [1486] "M00967_43_000000000-A3JHG_1_2102_25051_20557"
#> [1487] "M00967_43_000000000-A3JHG_1_2102_25569_21727"
#> [1488] "M00967_43_000000000-A3JHG_1_2102_25992_9452" 
#> [1489] "M00967_43_000000000-A3JHG_1_2102_26214_23057"
#> [1490] "M00967_43_000000000-A3JHG_1_2102_26372_21137"
#> [1491] "M00967_43_000000000-A3JHG_1_2102_26498_14617"
#> [1492] "M00967_43_000000000-A3JHG_1_2102_26655_15452"
#> [1493] "M00967_43_000000000-A3JHG_1_2102_27520_15637"
#> [1494] "M00967_43_000000000-A3JHG_1_2102_27720_18606"
#> [1495] "M00967_43_000000000-A3JHG_1_2102_27795_15713"
#> [1496] "M00967_43_000000000-A3JHG_1_2102_28446_16991"
#> [1497] "M00967_43_000000000-A3JHG_1_2102_5173_21862" 
#> [1498] "M00967_43_000000000-A3JHG_1_2102_5706_8003"  
#> [1499] "M00967_43_000000000-A3JHG_1_2102_5852_19432" 
#> [1500] "M00967_43_000000000-A3JHG_1_2102_5906_19319" 
#> [1501] "M00967_43_000000000-A3JHG_1_2102_6705_15205" 
#> [1502] "M00967_43_000000000-A3JHG_1_2102_7041_13746" 
#> [1503] "M00967_43_000000000-A3JHG_1_2102_7671_11355" 
#> [1504] "M00967_43_000000000-A3JHG_1_2102_7679_12922" 
#> [1505] "M00967_43_000000000-A3JHG_1_2102_7871_10547" 
#> [1506] "M00967_43_000000000-A3JHG_1_2102_8178_6298"  
#> [1507] "M00967_43_000000000-A3JHG_1_2102_8408_13436" 
#> [1508] "M00967_43_000000000-A3JHG_1_2102_8893_14291" 
#> [1509] "M00967_43_000000000-A3JHG_1_2102_9019_15362" 
#> [1510] "M00967_43_000000000-A3JHG_1_2102_9710_20744" 
#> [1511] "M00967_43_000000000-A3JHG_1_2103_10155_24625"
#> [1512] "M00967_43_000000000-A3JHG_1_2103_10530_5261" 
#> [1513] "M00967_43_000000000-A3JHG_1_2103_10569_14341"
#> [1514] "M00967_43_000000000-A3JHG_1_2103_10962_17123"
#> [1515] "M00967_43_000000000-A3JHG_1_2103_11158_26523"
#> [1516] "M00967_43_000000000-A3JHG_1_2103_11533_13455"
#> [1517] "M00967_43_000000000-A3JHG_1_2103_11591_8763" 
#> [1518] "M00967_43_000000000-A3JHG_1_2103_11635_15432"
#> [1519] "M00967_43_000000000-A3JHG_1_2103_11687_23017"
#> [1520] "M00967_43_000000000-A3JHG_1_2103_11783_23160"
#> [1521] "M00967_43_000000000-A3JHG_1_2103_11954_7925" 
#> [1522] "M00967_43_000000000-A3JHG_1_2103_12055_12246"
#> [1523] "M00967_43_000000000-A3JHG_1_2103_12193_17343"
#> [1524] "M00967_43_000000000-A3JHG_1_2103_12503_10753"
#> [1525] "M00967_43_000000000-A3JHG_1_2103_12547_24481"
#> [1526] "M00967_43_000000000-A3JHG_1_2103_12830_18492"
#> [1527] "M00967_43_000000000-A3JHG_1_2103_13084_5117" 
#> [1528] "M00967_43_000000000-A3JHG_1_2103_14180_15793"
#> [1529] "M00967_43_000000000-A3JHG_1_2103_14570_24653"
#> [1530] "M00967_43_000000000-A3JHG_1_2103_14978_21851"
#> [1531] "M00967_43_000000000-A3JHG_1_2103_15051_3502" 
#> [1532] "M00967_43_000000000-A3JHG_1_2103_15244_25817"
#> [1533] "M00967_43_000000000-A3JHG_1_2103_15275_25227"
#> [1534] "M00967_43_000000000-A3JHG_1_2103_15345_25846"
#> [1535] "M00967_43_000000000-A3JHG_1_2103_15482_15086"
#> [1536] "M00967_43_000000000-A3JHG_1_2103_15669_23137"
#> [1537] "M00967_43_000000000-A3JHG_1_2103_16331_20380"
#> [1538] "M00967_43_000000000-A3JHG_1_2103_16606_4023" 
#> [1539] "M00967_43_000000000-A3JHG_1_2103_16795_27764"
#> [1540] "M00967_43_000000000-A3JHG_1_2103_16979_24828"
#> [1541] "M00967_43_000000000-A3JHG_1_2103_17221_27256"
#> [1542] "M00967_43_000000000-A3JHG_1_2103_18492_19103"
#> [1543] "M00967_43_000000000-A3JHG_1_2103_18975_9767" 
#> [1544] "M00967_43_000000000-A3JHG_1_2103_19057_10430"
#> [1545] "M00967_43_000000000-A3JHG_1_2103_19125_26116"
#> [1546] "M00967_43_000000000-A3JHG_1_2103_19149_16619"
#> [1547] "M00967_43_000000000-A3JHG_1_2103_19578_6645" 
#> [1548] "M00967_43_000000000-A3JHG_1_2103_20113_16810"
#> [1549] "M00967_43_000000000-A3JHG_1_2103_20332_24524"
#> [1550] "M00967_43_000000000-A3JHG_1_2103_20753_25954"
#> [1551] "M00967_43_000000000-A3JHG_1_2103_21597_21834"
#> [1552] "M00967_43_000000000-A3JHG_1_2103_21648_16922"
#> [1553] "M00967_43_000000000-A3JHG_1_2103_21922_24386"
#> [1554] "M00967_43_000000000-A3JHG_1_2103_22029_15973"
#> [1555] "M00967_43_000000000-A3JHG_1_2103_22136_18992"
#> [1556] "M00967_43_000000000-A3JHG_1_2103_22380_17554"
#> [1557] "M00967_43_000000000-A3JHG_1_2103_22514_6355" 
#> [1558] "M00967_43_000000000-A3JHG_1_2103_23355_23413"
#> [1559] "M00967_43_000000000-A3JHG_1_2103_24284_21559"
#> [1560] "M00967_43_000000000-A3JHG_1_2103_24918_13396"
#> [1561] "M00967_43_000000000-A3JHG_1_2103_25452_6018" 
#> [1562] "M00967_43_000000000-A3JHG_1_2103_25809_24518"
#> [1563] "M00967_43_000000000-A3JHG_1_2103_25987_20682"
#> [1564] "M00967_43_000000000-A3JHG_1_2103_26632_14477"
#> [1565] "M00967_43_000000000-A3JHG_1_2103_26747_11121"
#> [1566] "M00967_43_000000000-A3JHG_1_2103_26975_11754"
#> [1567] "M00967_43_000000000-A3JHG_1_2103_27141_19101"
#> [1568] "M00967_43_000000000-A3JHG_1_2103_27613_20155"
#> [1569] "M00967_43_000000000-A3JHG_1_2103_27740_10743"
#> [1570] "M00967_43_000000000-A3JHG_1_2103_28700_17194"
#> [1571] "M00967_43_000000000-A3JHG_1_2103_29235_16481"
#> [1572] "M00967_43_000000000-A3JHG_1_2103_3688_9745"  
#> [1573] "M00967_43_000000000-A3JHG_1_2103_4125_17157" 
#> [1574] "M00967_43_000000000-A3JHG_1_2103_5527_9613"  
#> [1575] "M00967_43_000000000-A3JHG_1_2103_5532_14852" 
#> [1576] "M00967_43_000000000-A3JHG_1_2103_6243_11429" 
#> [1577] "M00967_43_000000000-A3JHG_1_2103_6640_14822" 
#> [1578] "M00967_43_000000000-A3JHG_1_2103_6845_19452" 
#> [1579] "M00967_43_000000000-A3JHG_1_2103_7701_13077" 
#> [1580] "M00967_43_000000000-A3JHG_1_2103_7735_17880" 
#> [1581] "M00967_43_000000000-A3JHG_1_2103_7853_22683" 
#> [1582] "M00967_43_000000000-A3JHG_1_2103_7856_11382" 
#> [1583] "M00967_43_000000000-A3JHG_1_2103_7887_6353"  
#> [1584] "M00967_43_000000000-A3JHG_1_2103_8111_16669" 
#> [1585] "M00967_43_000000000-A3JHG_1_2103_9090_14716" 
#> [1586] "M00967_43_000000000-A3JHG_1_2103_9222_10640" 
#> [1587] "M00967_43_000000000-A3JHG_1_2103_9278_23216" 
#> [1588] "M00967_43_000000000-A3JHG_1_2103_9542_23894" 
#> [1589] "M00967_43_000000000-A3JHG_1_2103_9558_19765" 
#> [1590] "M00967_43_000000000-A3JHG_1_2104_10372_21530"
#> [1591] "M00967_43_000000000-A3JHG_1_2104_10436_21721"
#> [1592] "M00967_43_000000000-A3JHG_1_2104_11018_16527"
#> [1593] "M00967_43_000000000-A3JHG_1_2104_11063_11086"
#> [1594] "M00967_43_000000000-A3JHG_1_2104_11759_16153"
#> [1595] "M00967_43_000000000-A3JHG_1_2104_11929_10165"
#> [1596] "M00967_43_000000000-A3JHG_1_2104_12385_21989"
#> [1597] "M00967_43_000000000-A3JHG_1_2104_12407_10383"
#> [1598] "M00967_43_000000000-A3JHG_1_2104_12505_3424" 
#> [1599] "M00967_43_000000000-A3JHG_1_2104_12560_11750"
#> [1600] "M00967_43_000000000-A3JHG_1_2104_13679_3191" 
#> [1601] "M00967_43_000000000-A3JHG_1_2104_13757_8649" 
#> [1602] "M00967_43_000000000-A3JHG_1_2104_14075_7093" 
#> [1603] "M00967_43_000000000-A3JHG_1_2104_14359_14505"
#> [1604] "M00967_43_000000000-A3JHG_1_2104_14573_20188"
#> [1605] "M00967_43_000000000-A3JHG_1_2104_14582_18801"
#> [1606] "M00967_43_000000000-A3JHG_1_2104_14864_25191"
#> [1607] "M00967_43_000000000-A3JHG_1_2104_15263_18242"
#> [1608] "M00967_43_000000000-A3JHG_1_2104_15288_13724"
#> [1609] "M00967_43_000000000-A3JHG_1_2104_15567_14710"
#> [1610] "M00967_43_000000000-A3JHG_1_2104_15684_7042" 
#> [1611] "M00967_43_000000000-A3JHG_1_2104_16003_12063"
#> [1612] "M00967_43_000000000-A3JHG_1_2104_17227_9135" 
#> [1613] "M00967_43_000000000-A3JHG_1_2104_17426_10838"
#> [1614] "M00967_43_000000000-A3JHG_1_2104_17794_25042"
#> [1615] "M00967_43_000000000-A3JHG_1_2104_17870_18881"
#> [1616] "M00967_43_000000000-A3JHG_1_2104_17991_24804"
#> [1617] "M00967_43_000000000-A3JHG_1_2104_18172_10876"
#> [1618] "M00967_43_000000000-A3JHG_1_2104_18217_19581"
#> [1619] "M00967_43_000000000-A3JHG_1_2104_18261_8446" 
#> [1620] "M00967_43_000000000-A3JHG_1_2104_18459_1792" 
#> [1621] "M00967_43_000000000-A3JHG_1_2104_18525_11294"
#> [1622] "M00967_43_000000000-A3JHG_1_2104_18899_8185" 
#> [1623] "M00967_43_000000000-A3JHG_1_2104_19264_18585"
#> [1624] "M00967_43_000000000-A3JHG_1_2104_19276_2046" 
#> [1625] "M00967_43_000000000-A3JHG_1_2104_19435_2222" 
#> [1626] "M00967_43_000000000-A3JHG_1_2104_19605_13841"
#> [1627] "M00967_43_000000000-A3JHG_1_2104_20064_5598" 
#> [1628] "M00967_43_000000000-A3JHG_1_2104_20074_25538"
#> [1629] "M00967_43_000000000-A3JHG_1_2104_20685_7431" 
#> [1630] "M00967_43_000000000-A3JHG_1_2104_21039_10813"
#> [1631] "M00967_43_000000000-A3JHG_1_2104_21537_12875"
#> [1632] "M00967_43_000000000-A3JHG_1_2104_21630_21997"
#> [1633] "M00967_43_000000000-A3JHG_1_2104_22154_4269" 
#> [1634] "M00967_43_000000000-A3JHG_1_2104_22645_18534"
#> [1635] "M00967_43_000000000-A3JHG_1_2104_23034_13380"
#> [1636] "M00967_43_000000000-A3JHG_1_2104_23298_12028"
#> [1637] "M00967_43_000000000-A3JHG_1_2104_23590_17131"
#> [1638] "M00967_43_000000000-A3JHG_1_2104_24158_9277" 
#> [1639] "M00967_43_000000000-A3JHG_1_2104_24218_17682"
#> [1640] "M00967_43_000000000-A3JHG_1_2104_24247_11022"
#> [1641] "M00967_43_000000000-A3JHG_1_2104_25226_10311"
#> [1642] "M00967_43_000000000-A3JHG_1_2104_25418_11257"
#> [1643] "M00967_43_000000000-A3JHG_1_2104_25610_6877" 
#> [1644] "M00967_43_000000000-A3JHG_1_2104_25692_5505" 
#> [1645] "M00967_43_000000000-A3JHG_1_2104_25789_8271" 
#> [1646] "M00967_43_000000000-A3JHG_1_2104_26095_18712"
#> [1647] "M00967_43_000000000-A3JHG_1_2104_26300_11226"
#> [1648] "M00967_43_000000000-A3JHG_1_2104_26311_10309"
#> [1649] "M00967_43_000000000-A3JHG_1_2104_26311_14104"
#> [1650] "M00967_43_000000000-A3JHG_1_2104_27131_15083"
#> [1651] "M00967_43_000000000-A3JHG_1_2104_27637_18457"
#> [1652] "M00967_43_000000000-A3JHG_1_2104_4696_10935" 
#> [1653] "M00967_43_000000000-A3JHG_1_2104_4792_21056" 
#> [1654] "M00967_43_000000000-A3JHG_1_2104_6580_11453" 
#> [1655] "M00967_43_000000000-A3JHG_1_2104_7819_20401" 
#> [1656] "M00967_43_000000000-A3JHG_1_2104_7953_16691" 
#> [1657] "M00967_43_000000000-A3JHG_1_2104_8006_7424"  
#> [1658] "M00967_43_000000000-A3JHG_1_2104_8408_16485" 
#> [1659] "M00967_43_000000000-A3JHG_1_2104_8663_16101" 
#> [1660] "M00967_43_000000000-A3JHG_1_2104_9368_14850" 
#> [1661] "M00967_43_000000000-A3JHG_1_2104_9824_7090"  
#> [1662] "M00967_43_000000000-A3JHG_1_2104_9943_6220"  
#> [1663] "M00967_43_000000000-A3JHG_1_2105_10122_15632"
#> [1664] "M00967_43_000000000-A3JHG_1_2105_10589_14796"
#> [1665] "M00967_43_000000000-A3JHG_1_2105_10725_6078" 
#> [1666] "M00967_43_000000000-A3JHG_1_2105_10730_3162" 
#> [1667] "M00967_43_000000000-A3JHG_1_2105_11723_5649" 
#> [1668] "M00967_43_000000000-A3JHG_1_2105_12113_5426" 
#> [1669] "M00967_43_000000000-A3JHG_1_2105_12118_8163" 
#> [1670] "M00967_43_000000000-A3JHG_1_2105_12511_24138"
#> [1671] "M00967_43_000000000-A3JHG_1_2105_12878_11098"
#> [1672] "M00967_43_000000000-A3JHG_1_2105_13267_15722"
#> [1673] "M00967_43_000000000-A3JHG_1_2105_13456_25322"
#> [1674] "M00967_43_000000000-A3JHG_1_2105_14589_9498" 
#> [1675] "M00967_43_000000000-A3JHG_1_2105_15118_19690"
#> [1676] "M00967_43_000000000-A3JHG_1_2105_15139_15368"
#> [1677] "M00967_43_000000000-A3JHG_1_2105_15512_14920"
#> [1678] "M00967_43_000000000-A3JHG_1_2105_16439_13484"
#> [1679] "M00967_43_000000000-A3JHG_1_2105_16818_5464" 
#> [1680] "M00967_43_000000000-A3JHG_1_2105_17039_19920"
#> [1681] "M00967_43_000000000-A3JHG_1_2105_17259_2092" 
#> [1682] "M00967_43_000000000-A3JHG_1_2105_17320_11249"
#> [1683] "M00967_43_000000000-A3JHG_1_2105_17326_26211"
#> [1684] "M00967_43_000000000-A3JHG_1_2105_17957_13444"
#> [1685] "M00967_43_000000000-A3JHG_1_2105_18070_5688" 
#> [1686] "M00967_43_000000000-A3JHG_1_2105_18677_14481"
#> [1687] "M00967_43_000000000-A3JHG_1_2105_18700_20485"
#> [1688] "M00967_43_000000000-A3JHG_1_2105_18816_26393"
#> [1689] "M00967_43_000000000-A3JHG_1_2105_19809_18745"
#> [1690] "M00967_43_000000000-A3JHG_1_2105_20695_7710" 
#> [1691] "M00967_43_000000000-A3JHG_1_2105_21230_21464"
#> [1692] "M00967_43_000000000-A3JHG_1_2105_21571_16049"
#> [1693] "M00967_43_000000000-A3JHG_1_2105_21901_10432"
#> [1694] "M00967_43_000000000-A3JHG_1_2105_22135_19871"
#> [1695] "M00967_43_000000000-A3JHG_1_2105_22686_4913" 
#> [1696] "M00967_43_000000000-A3JHG_1_2105_22749_25576"
#> [1697] "M00967_43_000000000-A3JHG_1_2105_23092_9359" 
#> [1698] "M00967_43_000000000-A3JHG_1_2105_23104_8731" 
#> [1699] "M00967_43_000000000-A3JHG_1_2105_23960_21462"
#> [1700] "M00967_43_000000000-A3JHG_1_2105_24305_5134" 
#> [1701] "M00967_43_000000000-A3JHG_1_2105_24670_14227"
#> [1702] "M00967_43_000000000-A3JHG_1_2105_24712_15313"
#> [1703] "M00967_43_000000000-A3JHG_1_2105_24795_12844"
#> [1704] "M00967_43_000000000-A3JHG_1_2105_25408_18555"
#> [1705] "M00967_43_000000000-A3JHG_1_2105_25509_16535"
#> [1706] "M00967_43_000000000-A3JHG_1_2105_25687_16182"
#> [1707] "M00967_43_000000000-A3JHG_1_2105_25960_21171"
#> [1708] "M00967_43_000000000-A3JHG_1_2105_26015_12877"
#> [1709] "M00967_43_000000000-A3JHG_1_2105_26628_10360"
#> [1710] "M00967_43_000000000-A3JHG_1_2105_27014_19197"
#> [1711] "M00967_43_000000000-A3JHG_1_2105_27616_16116"
#> [1712] "M00967_43_000000000-A3JHG_1_2105_29107_15046"
#> [1713] "M00967_43_000000000-A3JHG_1_2105_3878_16016" 
#> [1714] "M00967_43_000000000-A3JHG_1_2105_6171_20541" 
#> [1715] "M00967_43_000000000-A3JHG_1_2105_6844_12008" 
#> [1716] "M00967_43_000000000-A3JHG_1_2105_7116_11773" 
#> [1717] "M00967_43_000000000-A3JHG_1_2105_7169_20746" 
#> [1718] "M00967_43_000000000-A3JHG_1_2105_7943_8043"  
#> [1719] "M00967_43_000000000-A3JHG_1_2105_7983_21245" 
#> [1720] "M00967_43_000000000-A3JHG_1_2105_8155_16465" 
#> [1721] "M00967_43_000000000-A3JHG_1_2105_9580_7709"  
#> [1722] "M00967_43_000000000-A3JHG_1_2106_10304_14471"
#> [1723] "M00967_43_000000000-A3JHG_1_2106_10410_19621"
#> [1724] "M00967_43_000000000-A3JHG_1_2106_10758_21336"
#> [1725] "M00967_43_000000000-A3JHG_1_2106_10776_6955" 
#> [1726] "M00967_43_000000000-A3JHG_1_2106_10832_13902"
#> [1727] "M00967_43_000000000-A3JHG_1_2106_11005_4817" 
#> [1728] "M00967_43_000000000-A3JHG_1_2106_11123_12964"
#> [1729] "M00967_43_000000000-A3JHG_1_2106_12253_21736"
#> [1730] "M00967_43_000000000-A3JHG_1_2106_12340_14438"
#> [1731] "M00967_43_000000000-A3JHG_1_2106_12440_13627"
#> [1732] "M00967_43_000000000-A3JHG_1_2106_12968_23126"
#> [1733] "M00967_43_000000000-A3JHG_1_2106_13686_15465"
#> [1734] "M00967_43_000000000-A3JHG_1_2106_14031_23991"
#> [1735] "M00967_43_000000000-A3JHG_1_2106_14094_18415"
#> [1736] "M00967_43_000000000-A3JHG_1_2106_14236_13388"
#> [1737] "M00967_43_000000000-A3JHG_1_2106_14305_11884"
#> [1738] "M00967_43_000000000-A3JHG_1_2106_14743_23547"
#> [1739] "M00967_43_000000000-A3JHG_1_2106_14859_14822"
#> [1740] "M00967_43_000000000-A3JHG_1_2106_15051_17101"
#> [1741] "M00967_43_000000000-A3JHG_1_2106_15186_22937"
#> [1742] "M00967_43_000000000-A3JHG_1_2106_15211_24821"
#> [1743] "M00967_43_000000000-A3JHG_1_2106_15292_3687" 
#> [1744] "M00967_43_000000000-A3JHG_1_2106_15738_17268"
#> [1745] "M00967_43_000000000-A3JHG_1_2106_16059_20423"
#> [1746] "M00967_43_000000000-A3JHG_1_2106_16184_19905"
#> [1747] "M00967_43_000000000-A3JHG_1_2106_16419_20442"
#> [1748] "M00967_43_000000000-A3JHG_1_2106_16743_9576" 
#> [1749] "M00967_43_000000000-A3JHG_1_2106_17102_4744" 
#> [1750] "M00967_43_000000000-A3JHG_1_2106_17117_12389"
#> [1751] "M00967_43_000000000-A3JHG_1_2106_17332_22169"
#> [1752] "M00967_43_000000000-A3JHG_1_2106_17516_3721" 
#> [1753] "M00967_43_000000000-A3JHG_1_2106_17640_11873"
#> [1754] "M00967_43_000000000-A3JHG_1_2106_17671_25995"
#> [1755] "M00967_43_000000000-A3JHG_1_2106_17741_19018"
#> [1756] "M00967_43_000000000-A3JHG_1_2106_17935_25987"
#> [1757] "M00967_43_000000000-A3JHG_1_2106_18053_21392"
#> [1758] "M00967_43_000000000-A3JHG_1_2106_18324_28644"
#> [1759] "M00967_43_000000000-A3JHG_1_2106_18414_12416"
#> [1760] "M00967_43_000000000-A3JHG_1_2106_18656_26472"
#> [1761] "M00967_43_000000000-A3JHG_1_2106_19274_24278"
#> [1762] "M00967_43_000000000-A3JHG_1_2106_19379_4847" 
#> [1763] "M00967_43_000000000-A3JHG_1_2106_19750_8249" 
#> [1764] "M00967_43_000000000-A3JHG_1_2106_19951_19593"
#> [1765] "M00967_43_000000000-A3JHG_1_2106_20224_25526"
#> [1766] "M00967_43_000000000-A3JHG_1_2106_20458_8348" 
#> [1767] "M00967_43_000000000-A3JHG_1_2106_20527_14577"
#> [1768] "M00967_43_000000000-A3JHG_1_2106_20871_13331"
#> [1769] "M00967_43_000000000-A3JHG_1_2106_20954_22445"
#> [1770] "M00967_43_000000000-A3JHG_1_2106_21285_10820"
#> [1771] "M00967_43_000000000-A3JHG_1_2106_21459_16802"
#> [1772] "M00967_43_000000000-A3JHG_1_2106_21633_11578"
#> [1773] "M00967_43_000000000-A3JHG_1_2106_22696_14527"
#> [1774] "M00967_43_000000000-A3JHG_1_2106_22746_5360" 
#> [1775] "M00967_43_000000000-A3JHG_1_2106_22868_26146"
#> [1776] "M00967_43_000000000-A3JHG_1_2106_23085_20010"
#> [1777] "M00967_43_000000000-A3JHG_1_2106_23199_7043" 
#> [1778] "M00967_43_000000000-A3JHG_1_2106_23285_7924" 
#> [1779] "M00967_43_000000000-A3JHG_1_2106_23373_7365" 
#> [1780] "M00967_43_000000000-A3JHG_1_2106_25081_10002"
#> [1781] "M00967_43_000000000-A3JHG_1_2106_25337_11987"
#> [1782] "M00967_43_000000000-A3JHG_1_2106_25519_17332"
#> [1783] "M00967_43_000000000-A3JHG_1_2106_25719_12422"
#> [1784] "M00967_43_000000000-A3JHG_1_2106_27387_19532"
#> [1785] "M00967_43_000000000-A3JHG_1_2106_27480_14501"
#> [1786] "M00967_43_000000000-A3JHG_1_2106_28000_12820"
#> [1787] "M00967_43_000000000-A3JHG_1_2106_28239_16178"
#> [1788] "M00967_43_000000000-A3JHG_1_2106_28405_10515"
#> [1789] "M00967_43_000000000-A3JHG_1_2106_3727_12917" 
#> [1790] "M00967_43_000000000-A3JHG_1_2106_5509_18056" 
#> [1791] "M00967_43_000000000-A3JHG_1_2106_5552_10582" 
#> [1792] "M00967_43_000000000-A3JHG_1_2106_6366_20022" 
#> [1793] "M00967_43_000000000-A3JHG_1_2106_6909_9337"  
#> [1794] "M00967_43_000000000-A3JHG_1_2106_7250_19443" 
#> [1795] "M00967_43_000000000-A3JHG_1_2106_7936_8291"  
#> [1796] "M00967_43_000000000-A3JHG_1_2106_7943_14658" 
#> [1797] "M00967_43_000000000-A3JHG_1_2106_8030_23327" 
#> [1798] "M00967_43_000000000-A3JHG_1_2106_8679_15164" 
#> [1799] "M00967_43_000000000-A3JHG_1_2106_8736_3855"  
#> [1800] "M00967_43_000000000-A3JHG_1_2106_9506_14302" 
#> [1801] "M00967_43_000000000-A3JHG_1_2106_9777_15013" 
#> [1802] "M00967_43_000000000-A3JHG_1_2106_9918_20510" 
#> [1803] "M00967_43_000000000-A3JHG_1_2107_10624_8294" 
#> [1804] "M00967_43_000000000-A3JHG_1_2107_10803_20222"
#> [1805] "M00967_43_000000000-A3JHG_1_2107_10874_20271"
#> [1806] "M00967_43_000000000-A3JHG_1_2107_11226_20980"
#> [1807] "M00967_43_000000000-A3JHG_1_2107_11282_7357" 
#> [1808] "M00967_43_000000000-A3JHG_1_2107_11302_8835" 
#> [1809] "M00967_43_000000000-A3JHG_1_2107_11601_20781"
#> [1810] "M00967_43_000000000-A3JHG_1_2107_12056_5129" 
#> [1811] "M00967_43_000000000-A3JHG_1_2107_12316_24495"
#> [1812] "M00967_43_000000000-A3JHG_1_2107_12427_15408"
#> [1813] "M00967_43_000000000-A3JHG_1_2107_12502_13322"
#> [1814] "M00967_43_000000000-A3JHG_1_2107_12720_21761"
#> [1815] "M00967_43_000000000-A3JHG_1_2107_13406_20132"
#> [1816] "M00967_43_000000000-A3JHG_1_2107_13720_13200"
#> [1817] "M00967_43_000000000-A3JHG_1_2107_13773_19305"
#> [1818] "M00967_43_000000000-A3JHG_1_2107_14096_15016"
#> [1819] "M00967_43_000000000-A3JHG_1_2107_14104_6500" 
#> [1820] "M00967_43_000000000-A3JHG_1_2107_14468_14971"
#> [1821] "M00967_43_000000000-A3JHG_1_2107_15433_11976"
#> [1822] "M00967_43_000000000-A3JHG_1_2107_15770_4019" 
#> [1823] "M00967_43_000000000-A3JHG_1_2107_16534_14769"
#> [1824] "M00967_43_000000000-A3JHG_1_2107_17002_3531" 
#> [1825] "M00967_43_000000000-A3JHG_1_2107_17498_8085" 
#> [1826] "M00967_43_000000000-A3JHG_1_2107_17502_25142"
#> [1827] "M00967_43_000000000-A3JHG_1_2107_18284_2758" 
#> [1828] "M00967_43_000000000-A3JHG_1_2107_18354_8608" 
#> [1829] "M00967_43_000000000-A3JHG_1_2107_18723_18535"
#> [1830] "M00967_43_000000000-A3JHG_1_2107_19127_1886" 
#> [1831] "M00967_43_000000000-A3JHG_1_2107_19509_25673"
#> [1832] "M00967_43_000000000-A3JHG_1_2107_19676_23496"
#> [1833] "M00967_43_000000000-A3JHG_1_2107_20170_24454"
#> [1834] "M00967_43_000000000-A3JHG_1_2107_20457_11366"
#> [1835] "M00967_43_000000000-A3JHG_1_2107_20718_23357"
#> [1836] "M00967_43_000000000-A3JHG_1_2107_21370_16376"
#> [1837] "M00967_43_000000000-A3JHG_1_2107_21426_22343"
#> [1838] "M00967_43_000000000-A3JHG_1_2107_21758_14508"
#> [1839] "M00967_43_000000000-A3JHG_1_2107_21893_17093"
#> [1840] "M00967_43_000000000-A3JHG_1_2107_22192_18692"
#> [1841] "M00967_43_000000000-A3JHG_1_2107_22934_10355"
#> [1842] "M00967_43_000000000-A3JHG_1_2107_22959_20262"
#> [1843] "M00967_43_000000000-A3JHG_1_2107_23336_11976"
#> [1844] "M00967_43_000000000-A3JHG_1_2107_23464_17422"
#> [1845] "M00967_43_000000000-A3JHG_1_2107_23641_20458"
#> [1846] "M00967_43_000000000-A3JHG_1_2107_23657_12369"
#> [1847] "M00967_43_000000000-A3JHG_1_2107_23675_17555"
#> [1848] "M00967_43_000000000-A3JHG_1_2107_23831_12272"
#> [1849] "M00967_43_000000000-A3JHG_1_2107_24787_25663"
#> [1850] "M00967_43_000000000-A3JHG_1_2107_25487_6354" 
#> [1851] "M00967_43_000000000-A3JHG_1_2107_25612_9862" 
#> [1852] "M00967_43_000000000-A3JHG_1_2107_25731_13155"
#> [1853] "M00967_43_000000000-A3JHG_1_2107_26555_14955"
#> [1854] "M00967_43_000000000-A3JHG_1_2107_26817_14749"
#> [1855] "M00967_43_000000000-A3JHG_1_2107_27195_14074"
#> [1856] "M00967_43_000000000-A3JHG_1_2107_27604_15299"
#> [1857] "M00967_43_000000000-A3JHG_1_2107_27657_17343"
#> [1858] "M00967_43_000000000-A3JHG_1_2107_27757_15813"
#> [1859] "M00967_43_000000000-A3JHG_1_2107_28168_19224"
#> [1860] "M00967_43_000000000-A3JHG_1_2107_4030_17110" 
#> [1861] "M00967_43_000000000-A3JHG_1_2107_4339_17322" 
#> [1862] "M00967_43_000000000-A3JHG_1_2107_5816_12823" 
#> [1863] "M00967_43_000000000-A3JHG_1_2107_6821_13209" 
#> [1864] "M00967_43_000000000-A3JHG_1_2107_6956_16009" 
#> [1865] "M00967_43_000000000-A3JHG_1_2107_7207_7404"  
#> [1866] "M00967_43_000000000-A3JHG_1_2107_7268_7252"  
#> [1867] "M00967_43_000000000-A3JHG_1_2107_7483_16182" 
#> [1868] "M00967_43_000000000-A3JHG_1_2107_7538_13991" 
#> [1869] "M00967_43_000000000-A3JHG_1_2107_8153_5482"  
#> [1870] "M00967_43_000000000-A3JHG_1_2107_8240_20277" 
#> [1871] "M00967_43_000000000-A3JHG_1_2107_9330_10454" 
#> [1872] "M00967_43_000000000-A3JHG_1_2108_10109_11627"
#> [1873] "M00967_43_000000000-A3JHG_1_2108_10122_16816"
#> [1874] "M00967_43_000000000-A3JHG_1_2108_10549_22228"
#> [1875] "M00967_43_000000000-A3JHG_1_2108_10602_19971"
#> [1876] "M00967_43_000000000-A3JHG_1_2108_11790_16027"
#> [1877] "M00967_43_000000000-A3JHG_1_2108_12449_21256"
#> [1878] "M00967_43_000000000-A3JHG_1_2108_12490_16505"
#> [1879] "M00967_43_000000000-A3JHG_1_2108_12567_15271"
#> [1880] "M00967_43_000000000-A3JHG_1_2108_12620_23717"
#> [1881] "M00967_43_000000000-A3JHG_1_2108_12744_24921"
#> [1882] "M00967_43_000000000-A3JHG_1_2108_12907_13326"
#> [1883] "M00967_43_000000000-A3JHG_1_2108_12975_25037"
#> [1884] "M00967_43_000000000-A3JHG_1_2108_13254_24352"
#> [1885] "M00967_43_000000000-A3JHG_1_2108_13445_6827" 
#> [1886] "M00967_43_000000000-A3JHG_1_2108_13526_22864"
#> [1887] "M00967_43_000000000-A3JHG_1_2108_13984_9505" 
#> [1888] "M00967_43_000000000-A3JHG_1_2108_14452_7789" 
#> [1889] "M00967_43_000000000-A3JHG_1_2108_14707_9807" 
#> [1890] "M00967_43_000000000-A3JHG_1_2108_14981_25247"
#> [1891] "M00967_43_000000000-A3JHG_1_2108_15232_28922"
#> [1892] "M00967_43_000000000-A3JHG_1_2108_15244_28739"
#> [1893] "M00967_43_000000000-A3JHG_1_2108_15661_13758"
#> [1894] "M00967_43_000000000-A3JHG_1_2108_16172_11937"
#> [1895] "M00967_43_000000000-A3JHG_1_2108_16551_3455" 
#> [1896] "M00967_43_000000000-A3JHG_1_2108_16667_9805" 
#> [1897] "M00967_43_000000000-A3JHG_1_2108_16914_2576" 
#> [1898] "M00967_43_000000000-A3JHG_1_2108_17213_14575"
#> [1899] "M00967_43_000000000-A3JHG_1_2108_17448_14122"
#> [1900] "M00967_43_000000000-A3JHG_1_2108_17613_16563"
#> [1901] "M00967_43_000000000-A3JHG_1_2108_17961_17141"
#> [1902] "M00967_43_000000000-A3JHG_1_2108_18297_25797"
#> [1903] "M00967_43_000000000-A3JHG_1_2108_18302_14612"
#> [1904] "M00967_43_000000000-A3JHG_1_2108_18358_2270" 
#> [1905] "M00967_43_000000000-A3JHG_1_2108_18366_21641"
#> [1906] "M00967_43_000000000-A3JHG_1_2108_18470_11561"
#> [1907] "M00967_43_000000000-A3JHG_1_2108_18545_16099"
#> [1908] "M00967_43_000000000-A3JHG_1_2108_18841_25194"
#> [1909] "M00967_43_000000000-A3JHG_1_2108_18874_21200"
#> [1910] "M00967_43_000000000-A3JHG_1_2108_19062_16604"
#> [1911] "M00967_43_000000000-A3JHG_1_2108_19309_25091"
#> [1912] "M00967_43_000000000-A3JHG_1_2108_19840_26409"
#> [1913] "M00967_43_000000000-A3JHG_1_2108_19943_26491"
#> [1914] "M00967_43_000000000-A3JHG_1_2108_19974_4752" 
#> [1915] "M00967_43_000000000-A3JHG_1_2108_20440_11777"
#> [1916] "M00967_43_000000000-A3JHG_1_2108_20601_20789"
#> [1917] "M00967_43_000000000-A3JHG_1_2108_20757_5103" 
#> [1918] "M00967_43_000000000-A3JHG_1_2108_21552_11279"
#> [1919] "M00967_43_000000000-A3JHG_1_2108_21817_19761"
#> [1920] "M00967_43_000000000-A3JHG_1_2108_21978_6235" 
#> [1921] "M00967_43_000000000-A3JHG_1_2108_22020_21434"
#> [1922] "M00967_43_000000000-A3JHG_1_2108_22087_14118"
#> [1923] "M00967_43_000000000-A3JHG_1_2108_22189_18381"
#> [1924] "M00967_43_000000000-A3JHG_1_2108_22213_11738"
#> [1925] "M00967_43_000000000-A3JHG_1_2108_22507_11051"
#> [1926] "M00967_43_000000000-A3JHG_1_2108_22647_16398"
#> [1927] "M00967_43_000000000-A3JHG_1_2108_22703_14521"
#> [1928] "M00967_43_000000000-A3JHG_1_2108_23196_9131" 
#> [1929] "M00967_43_000000000-A3JHG_1_2108_23299_14833"
#> [1930] "M00967_43_000000000-A3JHG_1_2108_23551_13901"
#> [1931] "M00967_43_000000000-A3JHG_1_2108_23615_21129"
#> [1932] "M00967_43_000000000-A3JHG_1_2108_24082_19618"
#> [1933] "M00967_43_000000000-A3JHG_1_2108_24812_22312"
#> [1934] "M00967_43_000000000-A3JHG_1_2108_25681_15613"
#> [1935] "M00967_43_000000000-A3JHG_1_2108_26358_22652"
#> [1936] "M00967_43_000000000-A3JHG_1_2108_26742_16490"
#> [1937] "M00967_43_000000000-A3JHG_1_2108_27303_17375"
#> [1938] "M00967_43_000000000-A3JHG_1_2108_27497_14930"
#> [1939] "M00967_43_000000000-A3JHG_1_2108_27613_12769"
#> [1940] "M00967_43_000000000-A3JHG_1_2108_27773_22176"
#> [1941] "M00967_43_000000000-A3JHG_1_2108_27865_15182"
#> [1942] "M00967_43_000000000-A3JHG_1_2108_27871_9754" 
#> [1943] "M00967_43_000000000-A3JHG_1_2108_27907_12796"
#> [1944] "M00967_43_000000000-A3JHG_1_2108_28318_19174"
#> [1945] "M00967_43_000000000-A3JHG_1_2108_29278_13914"
#> [1946] "M00967_43_000000000-A3JHG_1_2108_4168_17142" 
#> [1947] "M00967_43_000000000-A3JHG_1_2108_4691_12067" 
#> [1948] "M00967_43_000000000-A3JHG_1_2108_5077_12621" 
#> [1949] "M00967_43_000000000-A3JHG_1_2108_5881_15408" 
#> [1950] "M00967_43_000000000-A3JHG_1_2108_5933_20471" 
#> [1951] "M00967_43_000000000-A3JHG_1_2108_6233_13205" 
#> [1952] "M00967_43_000000000-A3JHG_1_2108_6318_9727"  
#> [1953] "M00967_43_000000000-A3JHG_1_2108_7561_16231" 
#> [1954] "M00967_43_000000000-A3JHG_1_2108_7624_17316" 
#> [1955] "M00967_43_000000000-A3JHG_1_2108_7645_19778" 
#> [1956] "M00967_43_000000000-A3JHG_1_2108_7674_17990" 
#> [1957] "M00967_43_000000000-A3JHG_1_2108_7751_22068" 
#> [1958] "M00967_43_000000000-A3JHG_1_2108_7953_13304" 
#> [1959] "M00967_43_000000000-A3JHG_1_2108_8189_19210" 
#> [1960] "M00967_43_000000000-A3JHG_1_2108_8220_23281" 
#> [1961] "M00967_43_000000000-A3JHG_1_2108_9350_12884" 
#> [1962] "M00967_43_000000000-A3JHG_1_2109_10394_17608"
#> [1963] "M00967_43_000000000-A3JHG_1_2109_10445_7921" 
#> [1964] "M00967_43_000000000-A3JHG_1_2109_10631_18447"
#> [1965] "M00967_43_000000000-A3JHG_1_2109_10641_20956"
#> [1966] "M00967_43_000000000-A3JHG_1_2109_10681_7135" 
#> [1967] "M00967_43_000000000-A3JHG_1_2109_11177_4507" 
#> [1968] "M00967_43_000000000-A3JHG_1_2109_11269_19733"
#> [1969] "M00967_43_000000000-A3JHG_1_2109_11692_22443"
#> [1970] "M00967_43_000000000-A3JHG_1_2109_11755_2607" 
#> [1971] "M00967_43_000000000-A3JHG_1_2109_12163_4647" 
#> [1972] "M00967_43_000000000-A3JHG_1_2109_13011_10787"
#> [1973] "M00967_43_000000000-A3JHG_1_2109_13206_18588"
#> [1974] "M00967_43_000000000-A3JHG_1_2109_13274_20740"
#> [1975] "M00967_43_000000000-A3JHG_1_2109_13351_28335"
#> [1976] "M00967_43_000000000-A3JHG_1_2109_14117_4744" 
#> [1977] "M00967_43_000000000-A3JHG_1_2109_14322_20905"
#> [1978] "M00967_43_000000000-A3JHG_1_2109_14413_2449" 
#> [1979] "M00967_43_000000000-A3JHG_1_2109_14530_23810"
#> [1980] "M00967_43_000000000-A3JHG_1_2109_14780_27058"
#> [1981] "M00967_43_000000000-A3JHG_1_2109_14974_11618"
#> [1982] "M00967_43_000000000-A3JHG_1_2109_15065_8580" 
#> [1983] "M00967_43_000000000-A3JHG_1_2109_15158_9571" 
#> [1984] "M00967_43_000000000-A3JHG_1_2109_15602_23047"
#> [1985] "M00967_43_000000000-A3JHG_1_2109_15784_10393"
#> [1986] "M00967_43_000000000-A3JHG_1_2109_16258_1585" 
#> [1987] "M00967_43_000000000-A3JHG_1_2109_16433_12706"
#> [1988] "M00967_43_000000000-A3JHG_1_2109_16593_20511"
#> [1989] "M00967_43_000000000-A3JHG_1_2109_16956_10390"
#> [1990] "M00967_43_000000000-A3JHG_1_2109_17316_4550" 
#> [1991] "M00967_43_000000000-A3JHG_1_2109_17324_3447" 
#> [1992] "M00967_43_000000000-A3JHG_1_2109_17345_6668" 
#> [1993] "M00967_43_000000000-A3JHG_1_2109_17577_12359"
#> [1994] "M00967_43_000000000-A3JHG_1_2109_18108_6799" 
#> [1995] "M00967_43_000000000-A3JHG_1_2109_18300_12673"
#> [1996] "M00967_43_000000000-A3JHG_1_2109_18443_5251" 
#> [1997] "M00967_43_000000000-A3JHG_1_2109_19279_7836" 
#> [1998] "M00967_43_000000000-A3JHG_1_2109_19519_2320" 
#> [1999] "M00967_43_000000000-A3JHG_1_2109_19574_4320" 
#> [2000] "M00967_43_000000000-A3JHG_1_2109_19731_20280"
#> [2001] "M00967_43_000000000-A3JHG_1_2109_19808_24467"
#> [2002] "M00967_43_000000000-A3JHG_1_2109_19976_22044"
#> [2003] "M00967_43_000000000-A3JHG_1_2109_20085_20830"
#> [2004] "M00967_43_000000000-A3JHG_1_2109_20096_14197"
#> [2005] "M00967_43_000000000-A3JHG_1_2109_20252_6302" 
#> [2006] "M00967_43_000000000-A3JHG_1_2109_20493_16450"
#> [2007] "M00967_43_000000000-A3JHG_1_2109_20497_25999"
#> [2008] "M00967_43_000000000-A3JHG_1_2109_20708_17882"
#> [2009] "M00967_43_000000000-A3JHG_1_2109_21370_14363"
#> [2010] "M00967_43_000000000-A3JHG_1_2109_21549_9949" 
#> [2011] "M00967_43_000000000-A3JHG_1_2109_21713_13835"
#> [2012] "M00967_43_000000000-A3JHG_1_2109_22245_24849"
#> [2013] "M00967_43_000000000-A3JHG_1_2109_23117_11899"
#> [2014] "M00967_43_000000000-A3JHG_1_2109_23522_9652" 
#> [2015] "M00967_43_000000000-A3JHG_1_2109_24741_14194"
#> [2016] "M00967_43_000000000-A3JHG_1_2109_24829_22616"
#> [2017] "M00967_43_000000000-A3JHG_1_2109_24964_16600"
#> [2018] "M00967_43_000000000-A3JHG_1_2109_25048_7479" 
#> [2019] "M00967_43_000000000-A3JHG_1_2109_25056_22893"
#> [2020] "M00967_43_000000000-A3JHG_1_2109_25240_9047" 
#> [2021] "M00967_43_000000000-A3JHG_1_2109_25242_9758" 
#> [2022] "M00967_43_000000000-A3JHG_1_2109_25409_11127"
#> [2023] "M00967_43_000000000-A3JHG_1_2109_25546_10365"
#> [2024] "M00967_43_000000000-A3JHG_1_2109_25599_13567"
#> [2025] "M00967_43_000000000-A3JHG_1_2109_25923_16574"
#> [2026] "M00967_43_000000000-A3JHG_1_2109_27058_11840"
#> [2027] "M00967_43_000000000-A3JHG_1_2109_27247_19486"
#> [2028] "M00967_43_000000000-A3JHG_1_2109_27297_13184"
#> [2029] "M00967_43_000000000-A3JHG_1_2109_27908_9851" 
#> [2030] "M00967_43_000000000-A3JHG_1_2109_28564_10017"
#> [2031] "M00967_43_000000000-A3JHG_1_2109_29414_16449"
#> [2032] "M00967_43_000000000-A3JHG_1_2109_3291_16466" 
#> [2033] "M00967_43_000000000-A3JHG_1_2109_6453_7009"  
#> [2034] "M00967_43_000000000-A3JHG_1_2109_6640_9328"  
#> [2035] "M00967_43_000000000-A3JHG_1_2109_6753_10007" 
#> [2036] "M00967_43_000000000-A3JHG_1_2109_6882_14781" 
#> [2037] "M00967_43_000000000-A3JHG_1_2109_7315_17427" 
#> [2038] "M00967_43_000000000-A3JHG_1_2109_7813_4701"  
#> [2039] "M00967_43_000000000-A3JHG_1_2109_8100_6611"  
#> [2040] "M00967_43_000000000-A3JHG_1_2109_8160_14127" 
#> [2041] "M00967_43_000000000-A3JHG_1_2109_9362_22992" 
#> [2042] "M00967_43_000000000-A3JHG_1_2109_9416_24199" 
#> [2043] "M00967_43_000000000-A3JHG_1_2109_9657_11455" 
#> [2044] "M00967_43_000000000-A3JHG_1_2109_9943_10858" 
#> [2045] "M00967_43_000000000-A3JHG_1_2110_10010_14992"
#> [2046] "M00967_43_000000000-A3JHG_1_2110_10107_20150"
#> [2047] "M00967_43_000000000-A3JHG_1_2110_10506_11014"
#> [2048] "M00967_43_000000000-A3JHG_1_2110_10606_15947"
#> [2049] "M00967_43_000000000-A3JHG_1_2110_11195_7814" 
#> [2050] "M00967_43_000000000-A3JHG_1_2110_11463_19602"
#> [2051] "M00967_43_000000000-A3JHG_1_2110_11592_18531"
#> [2052] "M00967_43_000000000-A3JHG_1_2110_11932_12993"
#> [2053] "M00967_43_000000000-A3JHG_1_2110_13298_25616"
#> [2054] "M00967_43_000000000-A3JHG_1_2110_13356_9691" 
#> [2055] "M00967_43_000000000-A3JHG_1_2110_13599_17305"
#> [2056] "M00967_43_000000000-A3JHG_1_2110_14094_18100"
#> [2057] "M00967_43_000000000-A3JHG_1_2110_14767_15469"
#> [2058] "M00967_43_000000000-A3JHG_1_2110_15016_5261" 
#> [2059] "M00967_43_000000000-A3JHG_1_2110_15516_28587"
#> [2060] "M00967_43_000000000-A3JHG_1_2110_16087_14724"
#> [2061] "M00967_43_000000000-A3JHG_1_2110_16282_2095" 
#> [2062] "M00967_43_000000000-A3JHG_1_2110_16325_18649"
#> [2063] "M00967_43_000000000-A3JHG_1_2110_16565_10623"
#> [2064] "M00967_43_000000000-A3JHG_1_2110_16698_16415"
#> [2065] "M00967_43_000000000-A3JHG_1_2110_16815_14755"
#> [2066] "M00967_43_000000000-A3JHG_1_2110_16975_4772" 
#> [2067] "M00967_43_000000000-A3JHG_1_2110_17263_11534"
#> [2068] "M00967_43_000000000-A3JHG_1_2110_17584_23187"
#> [2069] "M00967_43_000000000-A3JHG_1_2110_17635_9975" 
#> [2070] "M00967_43_000000000-A3JHG_1_2110_17772_24108"
#> [2071] "M00967_43_000000000-A3JHG_1_2110_18060_5720" 
#> [2072] "M00967_43_000000000-A3JHG_1_2110_18835_13887"
#> [2073] "M00967_43_000000000-A3JHG_1_2110_18978_25958"
#> [2074] "M00967_43_000000000-A3JHG_1_2110_18983_13738"
#> [2075] "M00967_43_000000000-A3JHG_1_2110_19140_21811"
#> [2076] "M00967_43_000000000-A3JHG_1_2110_19314_26014"
#> [2077] "M00967_43_000000000-A3JHG_1_2110_19623_12472"
#> [2078] "M00967_43_000000000-A3JHG_1_2110_20054_3531" 
#> [2079] "M00967_43_000000000-A3JHG_1_2110_20937_6188" 
#> [2080] "M00967_43_000000000-A3JHG_1_2110_21191_16122"
#> [2081] "M00967_43_000000000-A3JHG_1_2110_21271_27854"
#> [2082] "M00967_43_000000000-A3JHG_1_2110_22006_25345"
#> [2083] "M00967_43_000000000-A3JHG_1_2110_22037_7673" 
#> [2084] "M00967_43_000000000-A3JHG_1_2110_22767_24294"
#> [2085] "M00967_43_000000000-A3JHG_1_2110_22861_3771" 
#> [2086] "M00967_43_000000000-A3JHG_1_2110_23067_17619"
#> [2087] "M00967_43_000000000-A3JHG_1_2110_23128_22332"
#> [2088] "M00967_43_000000000-A3JHG_1_2110_23478_21034"
#> [2089] "M00967_43_000000000-A3JHG_1_2110_23584_15948"
#> [2090] "M00967_43_000000000-A3JHG_1_2110_23721_15987"
#> [2091] "M00967_43_000000000-A3JHG_1_2110_23990_16750"
#> [2092] "M00967_43_000000000-A3JHG_1_2110_24102_23570"
#> [2093] "M00967_43_000000000-A3JHG_1_2110_24113_10963"
#> [2094] "M00967_43_000000000-A3JHG_1_2110_24147_19370"
#> [2095] "M00967_43_000000000-A3JHG_1_2110_24314_8984" 
#> [2096] "M00967_43_000000000-A3JHG_1_2110_24465_19572"
#> [2097] "M00967_43_000000000-A3JHG_1_2110_25176_12953"
#> [2098] "M00967_43_000000000-A3JHG_1_2110_25915_23958"
#> [2099] "M00967_43_000000000-A3JHG_1_2110_26207_14420"
#> [2100] "M00967_43_000000000-A3JHG_1_2110_26547_22871"
#> [2101] "M00967_43_000000000-A3JHG_1_2110_26743_23224"
#> [2102] "M00967_43_000000000-A3JHG_1_2110_27146_22098"
#> [2103] "M00967_43_000000000-A3JHG_1_2110_27183_12779"
#> [2104] "M00967_43_000000000-A3JHG_1_2110_27497_19046"
#> [2105] "M00967_43_000000000-A3JHG_1_2110_4035_19768" 
#> [2106] "M00967_43_000000000-A3JHG_1_2110_5187_18269" 
#> [2107] "M00967_43_000000000-A3JHG_1_2110_5497_7982"  
#> [2108] "M00967_43_000000000-A3JHG_1_2110_5588_24334" 
#> [2109] "M00967_43_000000000-A3JHG_1_2110_6116_13104" 
#> [2110] "M00967_43_000000000-A3JHG_1_2110_6133_17005" 
#> [2111] "M00967_43_000000000-A3JHG_1_2110_7067_15061" 
#> [2112] "M00967_43_000000000-A3JHG_1_2110_7457_8637"  
#> [2113] "M00967_43_000000000-A3JHG_1_2110_7469_16799" 
#> [2114] "M00967_43_000000000-A3JHG_1_2110_7822_24159" 
#> [2115] "M00967_43_000000000-A3JHG_1_2110_8189_7507"  
#> [2116] "M00967_43_000000000-A3JHG_1_2110_8652_13094" 
#> [2117] "M00967_43_000000000-A3JHG_1_2110_8800_5111"  
#> [2118] "M00967_43_000000000-A3JHG_1_2110_8933_4869"  
#> [2119] "M00967_43_000000000-A3JHG_1_2110_8957_20034" 
#> [2120] "M00967_43_000000000-A3JHG_1_2110_9069_16610" 
#> [2121] "M00967_43_000000000-A3JHG_1_2110_9142_17889" 
#> [2122] "M00967_43_000000000-A3JHG_1_2110_9717_11591" 
#> [2123] "M00967_43_000000000-A3JHG_1_2110_9991_8767"  
#> [2124] "M00967_43_000000000-A3JHG_1_2111_10081_26254"
#> [2125] "M00967_43_000000000-A3JHG_1_2111_10149_9721" 
#> [2126] "M00967_43_000000000-A3JHG_1_2111_10373_11883"
#> [2127] "M00967_43_000000000-A3JHG_1_2111_11460_16881"
#> [2128] "M00967_43_000000000-A3JHG_1_2111_12166_13729"
#> [2129] "M00967_43_000000000-A3JHG_1_2111_12192_22382"
#> [2130] "M00967_43_000000000-A3JHG_1_2111_12205_7359" 
#> [2131] "M00967_43_000000000-A3JHG_1_2111_12365_24314"
#> [2132] "M00967_43_000000000-A3JHG_1_2111_12367_24198"
#> [2133] "M00967_43_000000000-A3JHG_1_2111_12498_13580"
#> [2134] "M00967_43_000000000-A3JHG_1_2111_12821_2190" 
#> [2135] "M00967_43_000000000-A3JHG_1_2111_13555_13333"
#> [2136] "M00967_43_000000000-A3JHG_1_2111_13579_10781"
#> [2137] "M00967_43_000000000-A3JHG_1_2111_13741_6972" 
#> [2138] "M00967_43_000000000-A3JHG_1_2111_13819_4746" 
#> [2139] "M00967_43_000000000-A3JHG_1_2111_14170_11158"
#> [2140] "M00967_43_000000000-A3JHG_1_2111_14452_16768"
#> [2141] "M00967_43_000000000-A3JHG_1_2111_14544_23278"
#> [2142] "M00967_43_000000000-A3JHG_1_2111_15523_6202" 
#> [2143] "M00967_43_000000000-A3JHG_1_2111_15748_5366" 
#> [2144] "M00967_43_000000000-A3JHG_1_2111_16058_17833"
#> [2145] "M00967_43_000000000-A3JHG_1_2111_16323_3522" 
#> [2146] "M00967_43_000000000-A3JHG_1_2111_16365_1730" 
#> [2147] "M00967_43_000000000-A3JHG_1_2111_16376_7553" 
#> [2148] "M00967_43_000000000-A3JHG_1_2111_16922_20621"
#> [2149] "M00967_43_000000000-A3JHG_1_2111_17074_12360"
#> [2150] "M00967_43_000000000-A3JHG_1_2111_17223_18856"
#> [2151] "M00967_43_000000000-A3JHG_1_2111_17362_10105"
#> [2152] "M00967_43_000000000-A3JHG_1_2111_17688_14654"
#> [2153] "M00967_43_000000000-A3JHG_1_2111_18339_15747"
#> [2154] "M00967_43_000000000-A3JHG_1_2111_18507_22576"
#> [2155] "M00967_43_000000000-A3JHG_1_2111_18550_10076"
#> [2156] "M00967_43_000000000-A3JHG_1_2111_18844_14068"
#> [2157] "M00967_43_000000000-A3JHG_1_2111_19248_7551" 
#> [2158] "M00967_43_000000000-A3JHG_1_2111_19511_20526"
#> [2159] "M00967_43_000000000-A3JHG_1_2111_19706_18823"
#> [2160] "M00967_43_000000000-A3JHG_1_2111_20415_6804" 
#> [2161] "M00967_43_000000000-A3JHG_1_2111_20930_15299"
#> [2162] "M00967_43_000000000-A3JHG_1_2111_21066_7506" 
#> [2163] "M00967_43_000000000-A3JHG_1_2111_21247_15845"
#> [2164] "M00967_43_000000000-A3JHG_1_2111_21529_21161"
#> [2165] "M00967_43_000000000-A3JHG_1_2111_22131_23386"
#> [2166] "M00967_43_000000000-A3JHG_1_2111_22278_7396" 
#> [2167] "M00967_43_000000000-A3JHG_1_2111_22417_22191"
#> [2168] "M00967_43_000000000-A3JHG_1_2111_23260_18439"
#> [2169] "M00967_43_000000000-A3JHG_1_2111_25105_17913"
#> [2170] "M00967_43_000000000-A3JHG_1_2111_25171_15931"
#> [2171] "M00967_43_000000000-A3JHG_1_2111_25352_17924"
#> [2172] "M00967_43_000000000-A3JHG_1_2111_25642_13609"
#> [2173] "M00967_43_000000000-A3JHG_1_2111_25701_15655"
#> [2174] "M00967_43_000000000-A3JHG_1_2111_26393_21929"
#> [2175] "M00967_43_000000000-A3JHG_1_2111_27603_17495"
#> [2176] "M00967_43_000000000-A3JHG_1_2111_28308_17388"
#> [2177] "M00967_43_000000000-A3JHG_1_2111_28638_18003"
#> [2178] "M00967_43_000000000-A3JHG_1_2111_28641_13458"
#> [2179] "M00967_43_000000000-A3JHG_1_2111_29227_17392"
#> [2180] "M00967_43_000000000-A3JHG_1_2111_4007_20870" 
#> [2181] "M00967_43_000000000-A3JHG_1_2111_4279_22940" 
#> [2182] "M00967_43_000000000-A3JHG_1_2111_4404_12788" 
#> [2183] "M00967_43_000000000-A3JHG_1_2111_4814_18702" 
#> [2184] "M00967_43_000000000-A3JHG_1_2111_4991_9664"  
#> [2185] "M00967_43_000000000-A3JHG_1_2111_5139_21514" 
#> [2186] "M00967_43_000000000-A3JHG_1_2111_5639_10651" 
#> [2187] "M00967_43_000000000-A3JHG_1_2111_5664_6671"  
#> [2188] "M00967_43_000000000-A3JHG_1_2111_5740_11365" 
#> [2189] "M00967_43_000000000-A3JHG_1_2111_6039_9306"  
#> [2190] "M00967_43_000000000-A3JHG_1_2111_6114_18144" 
#> [2191] "M00967_43_000000000-A3JHG_1_2111_6592_20714" 
#> [2192] "M00967_43_000000000-A3JHG_1_2111_6782_8678"  
#> [2193] "M00967_43_000000000-A3JHG_1_2111_6887_10490" 
#> [2194] "M00967_43_000000000-A3JHG_1_2111_6981_7811"  
#> [2195] "M00967_43_000000000-A3JHG_1_2111_7333_20307" 
#> [2196] "M00967_43_000000000-A3JHG_1_2111_7971_11807" 
#> [2197] "M00967_43_000000000-A3JHG_1_2111_8004_16360" 
#> [2198] "M00967_43_000000000-A3JHG_1_2111_8291_8111"  
#> [2199] "M00967_43_000000000-A3JHG_1_2111_8304_3695"  
#> [2200] "M00967_43_000000000-A3JHG_1_2111_8444_24070" 
#> [2201] "M00967_43_000000000-A3JHG_1_2111_9619_23737" 
#> [2202] "M00967_43_000000000-A3JHG_1_2112_10008_8529" 
#> [2203] "M00967_43_000000000-A3JHG_1_2112_10013_23090"
#> [2204] "M00967_43_000000000-A3JHG_1_2112_10029_17181"
#> [2205] "M00967_43_000000000-A3JHG_1_2112_10826_19542"
#> [2206] "M00967_43_000000000-A3JHG_1_2112_11313_13873"
#> [2207] "M00967_43_000000000-A3JHG_1_2112_11352_13241"
#> [2208] "M00967_43_000000000-A3JHG_1_2112_11494_8595" 
#> [2209] "M00967_43_000000000-A3JHG_1_2112_11736_17751"
#> [2210] "M00967_43_000000000-A3JHG_1_2112_11899_4014" 
#> [2211] "M00967_43_000000000-A3JHG_1_2112_12022_9082" 
#> [2212] "M00967_43_000000000-A3JHG_1_2112_12287_24914"
#> [2213] "M00967_43_000000000-A3JHG_1_2112_12596_22917"
#> [2214] "M00967_43_000000000-A3JHG_1_2112_12631_18383"
#> [2215] "M00967_43_000000000-A3JHG_1_2112_12749_8826" 
#> [2216] "M00967_43_000000000-A3JHG_1_2112_12800_22855"
#> [2217] "M00967_43_000000000-A3JHG_1_2112_13415_4791" 
#> [2218] "M00967_43_000000000-A3JHG_1_2112_13463_16561"
#> [2219] "M00967_43_000000000-A3JHG_1_2112_13608_21597"
#> [2220] "M00967_43_000000000-A3JHG_1_2112_14726_9028" 
#> [2221] "M00967_43_000000000-A3JHG_1_2112_14790_1829" 
#> [2222] "M00967_43_000000000-A3JHG_1_2112_14835_3016" 
#> [2223] "M00967_43_000000000-A3JHG_1_2112_15581_17973"
#> [2224] "M00967_43_000000000-A3JHG_1_2112_15899_17037"
#> [2225] "M00967_43_000000000-A3JHG_1_2112_16774_9147" 
#> [2226] "M00967_43_000000000-A3JHG_1_2112_16997_10935"
#> [2227] "M00967_43_000000000-A3JHG_1_2112_17191_16655"
#> [2228] "M00967_43_000000000-A3JHG_1_2112_17219_14317"
#> [2229] "M00967_43_000000000-A3JHG_1_2112_17281_21576"
#> [2230] "M00967_43_000000000-A3JHG_1_2112_18144_24421"
#> [2231] "M00967_43_000000000-A3JHG_1_2112_18311_25065"
#> [2232] "M00967_43_000000000-A3JHG_1_2112_18393_28178"
#> [2233] "M00967_43_000000000-A3JHG_1_2112_19379_3300" 
#> [2234] "M00967_43_000000000-A3JHG_1_2112_19543_24449"
#> [2235] "M00967_43_000000000-A3JHG_1_2112_19844_23079"
#> [2236] "M00967_43_000000000-A3JHG_1_2112_20105_17530"
#> [2237] "M00967_43_000000000-A3JHG_1_2112_20183_27184"
#> [2238] "M00967_43_000000000-A3JHG_1_2112_20644_19641"
#> [2239] "M00967_43_000000000-A3JHG_1_2112_20678_18094"
#> [2240] "M00967_43_000000000-A3JHG_1_2112_21004_8857" 
#> [2241] "M00967_43_000000000-A3JHG_1_2112_21019_6342" 
#> [2242] "M00967_43_000000000-A3JHG_1_2112_21038_9589" 
#> [2243] "M00967_43_000000000-A3JHG_1_2112_21752_23067"
#> [2244] "M00967_43_000000000-A3JHG_1_2112_22227_7717" 
#> [2245] "M00967_43_000000000-A3JHG_1_2112_22323_22944"
#> [2246] "M00967_43_000000000-A3JHG_1_2112_22596_18948"
#> [2247] "M00967_43_000000000-A3JHG_1_2112_22672_12392"
#> [2248] "M00967_43_000000000-A3JHG_1_2112_22673_23252"
#> [2249] "M00967_43_000000000-A3JHG_1_2112_22880_7920" 
#> [2250] "M00967_43_000000000-A3JHG_1_2112_22954_15571"
#> [2251] "M00967_43_000000000-A3JHG_1_2112_23259_13531"
#> [2252] "M00967_43_000000000-A3JHG_1_2112_23298_24981"
#> [2253] "M00967_43_000000000-A3JHG_1_2112_23618_22699"
#> [2254] "M00967_43_000000000-A3JHG_1_2112_24110_22115"
#> [2255] "M00967_43_000000000-A3JHG_1_2112_24260_7949" 
#> [2256] "M00967_43_000000000-A3JHG_1_2112_24286_16039"
#> [2257] "M00967_43_000000000-A3JHG_1_2112_24625_19641"
#> [2258] "M00967_43_000000000-A3JHG_1_2112_24725_10550"
#> [2259] "M00967_43_000000000-A3JHG_1_2112_24855_7054" 
#> [2260] "M00967_43_000000000-A3JHG_1_2112_25145_5016" 
#> [2261] "M00967_43_000000000-A3JHG_1_2112_25256_24071"
#> [2262] "M00967_43_000000000-A3JHG_1_2112_25758_22797"
#> [2263] "M00967_43_000000000-A3JHG_1_2112_25988_9682" 
#> [2264] "M00967_43_000000000-A3JHG_1_2112_26066_19401"
#> [2265] "M00967_43_000000000-A3JHG_1_2112_26207_7695" 
#> [2266] "M00967_43_000000000-A3JHG_1_2112_26291_17310"
#> [2267] "M00967_43_000000000-A3JHG_1_2112_26347_9344" 
#> [2268] "M00967_43_000000000-A3JHG_1_2112_26672_16660"
#> [2269] "M00967_43_000000000-A3JHG_1_2112_26773_14017"
#> [2270] "M00967_43_000000000-A3JHG_1_2112_26969_15672"
#> [2271] "M00967_43_000000000-A3JHG_1_2112_2702_12868" 
#> [2272] "M00967_43_000000000-A3JHG_1_2112_27658_13780"
#> [2273] "M00967_43_000000000-A3JHG_1_2112_28271_11256"
#> [2274] "M00967_43_000000000-A3JHG_1_2112_3905_12418" 
#> [2275] "M00967_43_000000000-A3JHG_1_2112_4873_16721" 
#> [2276] "M00967_43_000000000-A3JHG_1_2112_6716_23007" 
#> [2277] "M00967_43_000000000-A3JHG_1_2112_7693_18206" 
#> [2278] "M00967_43_000000000-A3JHG_1_2112_7837_7502"  
#> [2279] "M00967_43_000000000-A3JHG_1_2112_7963_5767"  
#> [2280] "M00967_43_000000000-A3JHG_1_2112_8019_22340" 
#> [2281] "M00967_43_000000000-A3JHG_1_2112_8202_24000" 
#> [2282] "M00967_43_000000000-A3JHG_1_2112_8250_17966" 
#> [2283] "M00967_43_000000000-A3JHG_1_2112_8454_7758"  
#> [2284] "M00967_43_000000000-A3JHG_1_2112_8686_5552"  
#> [2285] "M00967_43_000000000-A3JHG_1_2112_8777_8861"  
#> [2286] "M00967_43_000000000-A3JHG_1_2112_9093_4236"  
#> [2287] "M00967_43_000000000-A3JHG_1_2112_9712_10612" 
#> [2288] "M00967_43_000000000-A3JHG_1_2112_9811_9982"  
#> [2289] "M00967_43_000000000-A3JHG_1_2113_10242_4487" 
#> [2290] "M00967_43_000000000-A3JHG_1_2113_10370_2988" 
#> [2291] "M00967_43_000000000-A3JHG_1_2113_10949_9667" 
#> [2292] "M00967_43_000000000-A3JHG_1_2113_11757_15493"
#> [2293] "M00967_43_000000000-A3JHG_1_2113_11875_9361" 
#> [2294] "M00967_43_000000000-A3JHG_1_2113_11924_15441"
#> [2295] "M00967_43_000000000-A3JHG_1_2113_12032_3782" 
#> [2296] "M00967_43_000000000-A3JHG_1_2113_12406_27672"
#> [2297] "M00967_43_000000000-A3JHG_1_2113_12698_4375" 
#> [2298] "M00967_43_000000000-A3JHG_1_2113_13133_16013"
#> [2299] "M00967_43_000000000-A3JHG_1_2113_13275_27560"
#> [2300] "M00967_43_000000000-A3JHG_1_2113_13414_6546" 
#> [2301] "M00967_43_000000000-A3JHG_1_2113_13419_4336" 
#> [2302] "M00967_43_000000000-A3JHG_1_2113_13714_22583"
#> [2303] "M00967_43_000000000-A3JHG_1_2113_13744_1953" 
#> [2304] "M00967_43_000000000-A3JHG_1_2113_13747_4412" 
#> [2305] "M00967_43_000000000-A3JHG_1_2113_14487_9130" 
#> [2306] "M00967_43_000000000-A3JHG_1_2113_14931_21754"
#> [2307] "M00967_43_000000000-A3JHG_1_2113_15093_8390" 
#> [2308] "M00967_43_000000000-A3JHG_1_2113_15115_2610" 
#> [2309] "M00967_43_000000000-A3JHG_1_2113_15752_12601"
#> [2310] "M00967_43_000000000-A3JHG_1_2113_15802_4098" 
#> [2311] "M00967_43_000000000-A3JHG_1_2113_16554_11081"
#> [2312] "M00967_43_000000000-A3JHG_1_2113_16796_6663" 
#> [2313] "M00967_43_000000000-A3JHG_1_2113_16978_21215"
#> [2314] "M00967_43_000000000-A3JHG_1_2113_17952_9366" 
#> [2315] "M00967_43_000000000-A3JHG_1_2113_18085_11993"
#> [2316] "M00967_43_000000000-A3JHG_1_2113_18276_27021"
#> [2317] "M00967_43_000000000-A3JHG_1_2113_18397_7516" 
#> [2318] "M00967_43_000000000-A3JHG_1_2113_18844_7914" 
#> [2319] "M00967_43_000000000-A3JHG_1_2113_19864_14092"
#> [2320] "M00967_43_000000000-A3JHG_1_2113_20030_3857" 
#> [2321] "M00967_43_000000000-A3JHG_1_2113_20184_16464"
#> [2322] "M00967_43_000000000-A3JHG_1_2113_20319_11310"
#> [2323] "M00967_43_000000000-A3JHG_1_2113_20783_5301" 
#> [2324] "M00967_43_000000000-A3JHG_1_2113_20973_10900"
#> [2325] "M00967_43_000000000-A3JHG_1_2113_21028_10433"
#> [2326] "M00967_43_000000000-A3JHG_1_2113_21308_17042"
#> [2327] "M00967_43_000000000-A3JHG_1_2113_21352_15666"
#> [2328] "M00967_43_000000000-A3JHG_1_2113_21647_20109"
#> [2329] "M00967_43_000000000-A3JHG_1_2113_21662_22858"
#> [2330] "M00967_43_000000000-A3JHG_1_2113_21823_22689"
#> [2331] "M00967_43_000000000-A3JHG_1_2113_21899_25122"
#> [2332] "M00967_43_000000000-A3JHG_1_2113_22451_22884"
#> [2333] "M00967_43_000000000-A3JHG_1_2113_23168_9750" 
#> [2334] "M00967_43_000000000-A3JHG_1_2113_23466_13467"
#> [2335] "M00967_43_000000000-A3JHG_1_2113_23687_14507"
#> [2336] "M00967_43_000000000-A3JHG_1_2113_23848_13201"
#> [2337] "M00967_43_000000000-A3JHG_1_2113_24485_10653"
#> [2338] "M00967_43_000000000-A3JHG_1_2113_24996_25271"
#> [2339] "M00967_43_000000000-A3JHG_1_2113_25427_14638"
#> [2340] "M00967_43_000000000-A3JHG_1_2113_25465_9646" 
#> [2341] "M00967_43_000000000-A3JHG_1_2113_26675_15471"
#> [2342] "M00967_43_000000000-A3JHG_1_2113_27238_12024"
#> [2343] "M00967_43_000000000-A3JHG_1_2113_27299_8066" 
#> [2344] "M00967_43_000000000-A3JHG_1_2113_28097_13665"
#> [2345] "M00967_43_000000000-A3JHG_1_2113_29656_14155"
#> [2346] "M00967_43_000000000-A3JHG_1_2113_4057_19587" 
#> [2347] "M00967_43_000000000-A3JHG_1_2113_4525_9290"  
#> [2348] "M00967_43_000000000-A3JHG_1_2113_5333_12504" 
#> [2349] "M00967_43_000000000-A3JHG_1_2113_5393_13669" 
#> [2350] "M00967_43_000000000-A3JHG_1_2113_5453_21104" 
#> [2351] "M00967_43_000000000-A3JHG_1_2113_6545_8894"  
#> [2352] "M00967_43_000000000-A3JHG_1_2113_6721_20524" 
#> [2353] "M00967_43_000000000-A3JHG_1_2113_7021_18070" 
#> [2354] "M00967_43_000000000-A3JHG_1_2113_7083_22842" 
#> [2355] "M00967_43_000000000-A3JHG_1_2113_7165_8175"  
#> [2356] "M00967_43_000000000-A3JHG_1_2113_7343_26333" 
#> [2357] "M00967_43_000000000-A3JHG_1_2113_7430_24285" 
#> [2358] "M00967_43_000000000-A3JHG_1_2113_8012_12345" 
#> [2359] "M00967_43_000000000-A3JHG_1_2113_8502_12671" 
#> [2360] "M00967_43_000000000-A3JHG_1_2113_9492_17806" 
#> [2361] "M00967_43_000000000-A3JHG_1_2113_9583_13925" 
#> [2362] "M00967_43_000000000-A3JHG_1_2113_9708_9804"  
#> [2363] "M00967_43_000000000-A3JHG_1_2114_10056_13609"
#> [2364] "M00967_43_000000000-A3JHG_1_2114_10492_14599"
#> [2365] "M00967_43_000000000-A3JHG_1_2114_10524_26952"
#> [2366] "M00967_43_000000000-A3JHG_1_2114_11245_12909"
#> [2367] "M00967_43_000000000-A3JHG_1_2114_11892_8192" 
#> [2368] "M00967_43_000000000-A3JHG_1_2114_12364_4557" 
#> [2369] "M00967_43_000000000-A3JHG_1_2114_12423_13578"
#> [2370] "M00967_43_000000000-A3JHG_1_2114_12634_10967"
#> [2371] "M00967_43_000000000-A3JHG_1_2114_12741_23025"
#> [2372] "M00967_43_000000000-A3JHG_1_2114_13345_25809"
#> [2373] "M00967_43_000000000-A3JHG_1_2114_13884_15208"
#> [2374] "M00967_43_000000000-A3JHG_1_2114_14201_8661" 
#> [2375] "M00967_43_000000000-A3JHG_1_2114_14273_14388"
#> [2376] "M00967_43_000000000-A3JHG_1_2114_14906_9268" 
#> [2377] "M00967_43_000000000-A3JHG_1_2114_15073_13278"
#> [2378] "M00967_43_000000000-A3JHG_1_2114_15082_26647"
#> [2379] "M00967_43_000000000-A3JHG_1_2114_15143_1651" 
#> [2380] "M00967_43_000000000-A3JHG_1_2114_16137_10353"
#> [2381] "M00967_43_000000000-A3JHG_1_2114_16268_21350"
#> [2382] "M00967_43_000000000-A3JHG_1_2114_16615_9462" 
#> [2383] "M00967_43_000000000-A3JHG_1_2114_16946_3813" 
#> [2384] "M00967_43_000000000-A3JHG_1_2114_17062_26446"
#> [2385] "M00967_43_000000000-A3JHG_1_2114_17435_3730" 
#> [2386] "M00967_43_000000000-A3JHG_1_2114_17684_11765"
#> [2387] "M00967_43_000000000-A3JHG_1_2114_18085_25907"
#> [2388] "M00967_43_000000000-A3JHG_1_2114_18278_17470"
#> [2389] "M00967_43_000000000-A3JHG_1_2114_18507_11303"
#> [2390] "M00967_43_000000000-A3JHG_1_2114_18555_17761"
#> [2391] "M00967_43_000000000-A3JHG_1_2114_18715_18416"
#> [2392] "M00967_43_000000000-A3JHG_1_2114_18928_10051"
#> [2393] "M00967_43_000000000-A3JHG_1_2114_19251_6041" 
#> [2394] "M00967_43_000000000-A3JHG_1_2114_19661_25608"
#> [2395] "M00967_43_000000000-A3JHG_1_2114_20048_11075"
#> [2396] "M00967_43_000000000-A3JHG_1_2114_20048_7664" 
#> [2397] "M00967_43_000000000-A3JHG_1_2114_20216_27116"
#> [2398] "M00967_43_000000000-A3JHG_1_2114_20314_10028"
#> [2399] "M00967_43_000000000-A3JHG_1_2114_20868_11538"
#> [2400] "M00967_43_000000000-A3JHG_1_2114_20876_7767" 
#> [2401] "M00967_43_000000000-A3JHG_1_2114_20979_14273"
#> [2402] "M00967_43_000000000-A3JHG_1_2114_21674_11680"
#> [2403] "M00967_43_000000000-A3JHG_1_2114_22468_4333" 
#> [2404] "M00967_43_000000000-A3JHG_1_2114_22474_15335"
#> [2405] "M00967_43_000000000-A3JHG_1_2114_24274_20575"
#> [2406] "M00967_43_000000000-A3JHG_1_2114_24618_11376"
#> [2407] "M00967_43_000000000-A3JHG_1_2114_24635_10372"
#> [2408] "M00967_43_000000000-A3JHG_1_2114_24873_12639"
#> [2409] "M00967_43_000000000-A3JHG_1_2114_25121_9113" 
#> [2410] "M00967_43_000000000-A3JHG_1_2114_26629_16356"
#> [2411] "M00967_43_000000000-A3JHG_1_2114_27537_17929"
#> [2412] "M00967_43_000000000-A3JHG_1_2114_28476_17852"
#> [2413] "M00967_43_000000000-A3JHG_1_2114_4298_19710" 
#> [2414] "M00967_43_000000000-A3JHG_1_2114_4499_18914" 
#> [2415] "M00967_43_000000000-A3JHG_1_2114_4609_11919" 
#> [2416] "M00967_43_000000000-A3JHG_1_2114_4725_7023"  
#> [2417] "M00967_43_000000000-A3JHG_1_2114_4913_11992" 
#> [2418] "M00967_43_000000000-A3JHG_1_2114_5448_14619" 
#> [2419] "M00967_43_000000000-A3JHG_1_2114_5702_10845" 
#> [2420] "M00967_43_000000000-A3JHG_1_2114_6259_8307"  
#> [2421] "M00967_43_000000000-A3JHG_1_2114_6705_14284" 
#> [2422] "M00967_43_000000000-A3JHG_1_2114_7386_6808"  
#> [2423] "M00967_43_000000000-A3JHG_1_2114_7946_17808" 
#> [2424] "M00967_43_000000000-A3JHG_1_2114_8962_26729" 
#> [2425] "M00967_43_000000000-A3JHG_1_2114_9526_8096"  

# To get the names of the sequences present sample 'F3D0'
miseq$names(type = "sequence", samples = c("F3D0"))
#>   [1] "M00967_43_000000000-A3JHG_1_1101_10133_8460" 
#>   [2] "M00967_43_000000000-A3JHG_1_1101_10331_23332"
#>   [3] "M00967_43_000000000-A3JHG_1_1101_11035_15765"
#>   [4] "M00967_43_000000000-A3JHG_1_1101_14364_8401" 
#>   [5] "M00967_43_000000000-A3JHG_1_1101_15591_4696" 
#>   [6] "M00967_43_000000000-A3JHG_1_1101_18143_13375"
#>   [7] "M00967_43_000000000-A3JHG_1_1101_18278_3345" 
#>   [8] "M00967_43_000000000-A3JHG_1_1101_18346_24737"
#>   [9] "M00967_43_000000000-A3JHG_1_1101_18922_4934" 
#>  [10] "M00967_43_000000000-A3JHG_1_1101_19534_17052"
#>  [11] "M00967_43_000000000-A3JHG_1_1101_21616_8560" 
#>  [12] "M00967_43_000000000-A3JHG_1_1101_23238_24359"
#>  [13] "M00967_43_000000000-A3JHG_1_1101_6836_23417" 
#>  [14] "M00967_43_000000000-A3JHG_1_1101_8868_4602"  
#>  [15] "M00967_43_000000000-A3JHG_1_1101_9331_5806"  
#>  [16] "M00967_43_000000000-A3JHG_1_1101_9553_14094" 
#>  [17] "M00967_43_000000000-A3JHG_1_1101_9620_19745" 
#>  [18] "M00967_43_000000000-A3JHG_1_1102_11115_19075"
#>  [19] "M00967_43_000000000-A3JHG_1_1102_14905_6542" 
#>  [20] "M00967_43_000000000-A3JHG_1_1102_15058_11924"
#>  [21] "M00967_43_000000000-A3JHG_1_1102_16785_23687"
#>  [22] "M00967_43_000000000-A3JHG_1_1102_17763_13764"
#>  [23] "M00967_43_000000000-A3JHG_1_1102_2114_15227" 
#>  [24] "M00967_43_000000000-A3JHG_1_1103_13364_5496" 
#>  [25] "M00967_43_000000000-A3JHG_1_1103_13966_3813" 
#>  [26] "M00967_43_000000000-A3JHG_1_1103_14518_16099"
#>  [27] "M00967_43_000000000-A3JHG_1_1103_19870_21567"
#>  [28] "M00967_43_000000000-A3JHG_1_1103_21051_5371" 
#>  [29] "M00967_43_000000000-A3JHG_1_1103_22490_21890"
#>  [30] "M00967_43_000000000-A3JHG_1_1103_23102_6138" 
#>  [31] "M00967_43_000000000-A3JHG_1_1103_24429_8227" 
#>  [32] "M00967_43_000000000-A3JHG_1_1103_25888_23585"
#>  [33] "M00967_43_000000000-A3JHG_1_1103_3864_17599" 
#>  [34] "M00967_43_000000000-A3JHG_1_1103_5171_14027" 
#>  [35] "M00967_43_000000000-A3JHG_1_1103_5754_14689" 
#>  [36] "M00967_43_000000000-A3JHG_1_1103_7049_5487"  
#>  [37] "M00967_43_000000000-A3JHG_1_1103_8900_10891" 
#>  [38] "M00967_43_000000000-A3JHG_1_1103_9123_12021" 
#>  [39] "M00967_43_000000000-A3JHG_1_1104_10834_8009" 
#>  [40] "M00967_43_000000000-A3JHG_1_1104_14252_27883"
#>  [41] "M00967_43_000000000-A3JHG_1_1104_17583_11048"
#>  [42] "M00967_43_000000000-A3JHG_1_1104_19441_6353" 
#>  [43] "M00967_43_000000000-A3JHG_1_1104_25962_6708" 
#>  [44] "M00967_43_000000000-A3JHG_1_1105_11964_4686" 
#>  [45] "M00967_43_000000000-A3JHG_1_1105_12761_22353"
#>  [46] "M00967_43_000000000-A3JHG_1_1105_13657_20414"
#>  [47] "M00967_43_000000000-A3JHG_1_1105_14547_12843"
#>  [48] "M00967_43_000000000-A3JHG_1_1105_14835_21938"
#>  [49] "M00967_43_000000000-A3JHG_1_1105_18475_22342"
#>  [50] "M00967_43_000000000-A3JHG_1_1105_18566_14540"
#>  [51] "M00967_43_000000000-A3JHG_1_1105_22508_16792"
#>  [52] "M00967_43_000000000-A3JHG_1_1105_23025_12058"
#>  [53] "M00967_43_000000000-A3JHG_1_1105_23612_18072"
#>  [54] "M00967_43_000000000-A3JHG_1_1105_25421_15745"
#>  [55] "M00967_43_000000000-A3JHG_1_1105_25554_16982"
#>  [56] "M00967_43_000000000-A3JHG_1_1105_26683_22024"
#>  [57] "M00967_43_000000000-A3JHG_1_1105_5158_15329" 
#>  [58] "M00967_43_000000000-A3JHG_1_1105_9179_20196" 
#>  [59] "M00967_43_000000000-A3JHG_1_1106_11240_22282"
#>  [60] "M00967_43_000000000-A3JHG_1_1106_12231_13452"
#>  [61] "M00967_43_000000000-A3JHG_1_1106_20287_11044"
#>  [62] "M00967_43_000000000-A3JHG_1_1106_21375_27109"
#>  [63] "M00967_43_000000000-A3JHG_1_1106_26087_11342"
#>  [64] "M00967_43_000000000-A3JHG_1_1106_26888_10548"
#>  [65] "M00967_43_000000000-A3JHG_1_1106_9558_4409"  
#>  [66] "M00967_43_000000000-A3JHG_1_1107_10661_18652"
#>  [67] "M00967_43_000000000-A3JHG_1_1107_12904_20713"
#>  [68] "M00967_43_000000000-A3JHG_1_1107_13192_19802"
#>  [69] "M00967_43_000000000-A3JHG_1_1107_14548_18902"
#>  [70] "M00967_43_000000000-A3JHG_1_1107_15750_18592"
#>  [71] "M00967_43_000000000-A3JHG_1_1107_18594_19714"
#>  [72] "M00967_43_000000000-A3JHG_1_1107_18662_2044" 
#>  [73] "M00967_43_000000000-A3JHG_1_1107_19005_21895"
#>  [74] "M00967_43_000000000-A3JHG_1_1107_19088_18601"
#>  [75] "M00967_43_000000000-A3JHG_1_1107_19562_20695"
#>  [76] "M00967_43_000000000-A3JHG_1_1107_21038_8569" 
#>  [77] "M00967_43_000000000-A3JHG_1_1107_24854_21229"
#>  [78] "M00967_43_000000000-A3JHG_1_1107_6417_21653" 
#>  [79] "M00967_43_000000000-A3JHG_1_1107_7945_19118" 
#>  [80] "M00967_43_000000000-A3JHG_1_1108_11048_4561" 
#>  [81] "M00967_43_000000000-A3JHG_1_1108_12037_3053" 
#>  [82] "M00967_43_000000000-A3JHG_1_1108_15903_3841" 
#>  [83] "M00967_43_000000000-A3JHG_1_1108_18362_23908"
#>  [84] "M00967_43_000000000-A3JHG_1_1108_18996_9436" 
#>  [85] "M00967_43_000000000-A3JHG_1_1108_19065_7362" 
#>  [86] "M00967_43_000000000-A3JHG_1_1108_4844_10914" 
#>  [87] "M00967_43_000000000-A3JHG_1_1109_12679_21321"
#>  [88] "M00967_43_000000000-A3JHG_1_1109_13330_21597"
#>  [89] "M00967_43_000000000-A3JHG_1_1109_14359_12484"
#>  [90] "M00967_43_000000000-A3JHG_1_1109_15029_15883"
#>  [91] "M00967_43_000000000-A3JHG_1_1109_20401_24219"
#>  [92] "M00967_43_000000000-A3JHG_1_1109_20716_17379"
#>  [93] "M00967_43_000000000-A3JHG_1_1109_22380_20532"
#>  [94] "M00967_43_000000000-A3JHG_1_1109_24274_5733" 
#>  [95] "M00967_43_000000000-A3JHG_1_1109_26108_17578"
#>  [96] "M00967_43_000000000-A3JHG_1_1109_7971_24259" 
#>  [97] "M00967_43_000000000-A3JHG_1_1109_8014_5010"  
#>  [98] "M00967_43_000000000-A3JHG_1_1109_9629_6769"  
#>  [99] "M00967_43_000000000-A3JHG_1_1110_12295_19699"
#> [100] "M00967_43_000000000-A3JHG_1_1110_13033_21990"
#> [101] "M00967_43_000000000-A3JHG_1_1110_13719_21639"
#> [102] "M00967_43_000000000-A3JHG_1_1110_14443_28339"
#> [103] "M00967_43_000000000-A3JHG_1_1110_14717_26485"
#> [104] "M00967_43_000000000-A3JHG_1_1110_15399_7509" 
#> [105] "M00967_43_000000000-A3JHG_1_1110_17220_22884"
#> [106] "M00967_43_000000000-A3JHG_1_1110_18313_3175" 
#> [107] "M00967_43_000000000-A3JHG_1_1110_21542_3012" 
#> [108] "M00967_43_000000000-A3JHG_1_1110_5315_13833" 
#> [109] "M00967_43_000000000-A3JHG_1_1110_6952_6984"  
#> [110] "M00967_43_000000000-A3JHG_1_1111_11386_7041" 
#> [111] "M00967_43_000000000-A3JHG_1_1111_13795_8174" 
#> [112] "M00967_43_000000000-A3JHG_1_1111_17002_7999" 
#> [113] "M00967_43_000000000-A3JHG_1_1111_17819_12937"
#> [114] "M00967_43_000000000-A3JHG_1_1111_20588_10279"
#> [115] "M00967_43_000000000-A3JHG_1_1111_21761_7263" 
#> [116] "M00967_43_000000000-A3JHG_1_1111_22910_7741" 
#> [117] "M00967_43_000000000-A3JHG_1_1111_4893_9943"  
#> [118] "M00967_43_000000000-A3JHG_1_1111_7157_8555"  
#> [119] "M00967_43_000000000-A3JHG_1_1111_8373_10700" 
#> [120] "M00967_43_000000000-A3JHG_1_1111_8697_7063"  
#> [121] "M00967_43_000000000-A3JHG_1_1111_9294_7319"  
#> [122] "M00967_43_000000000-A3JHG_1_1111_9420_7291"  
#> [123] "M00967_43_000000000-A3JHG_1_1112_10024_21424"
#> [124] "M00967_43_000000000-A3JHG_1_1112_10194_20788"
#> [125] "M00967_43_000000000-A3JHG_1_1112_11192_11762"
#> [126] "M00967_43_000000000-A3JHG_1_1112_11343_21537"
#> [127] "M00967_43_000000000-A3JHG_1_1112_12715_8487" 
#> [128] "M00967_43_000000000-A3JHG_1_1112_12844_3769" 
#> [129] "M00967_43_000000000-A3JHG_1_1112_13142_18436"
#> [130] "M00967_43_000000000-A3JHG_1_1112_14780_8078" 
#> [131] "M00967_43_000000000-A3JHG_1_1112_17696_15281"
#> [132] "M00967_43_000000000-A3JHG_1_1112_18848_19266"
#> [133] "M00967_43_000000000-A3JHG_1_1112_19376_17765"
#> [134] "M00967_43_000000000-A3JHG_1_1112_20127_11122"
#> [135] "M00967_43_000000000-A3JHG_1_1112_20884_13058"
#> [136] "M00967_43_000000000-A3JHG_1_1112_24118_19062"
#> [137] "M00967_43_000000000-A3JHG_1_1112_24606_18511"
#> [138] "M00967_43_000000000-A3JHG_1_1112_25719_18946"
#> [139] "M00967_43_000000000-A3JHG_1_1112_6862_18037" 
#> [140] "M00967_43_000000000-A3JHG_1_1112_7075_7100"  
#> [141] "M00967_43_000000000-A3JHG_1_1112_8821_13863" 
#> [142] "M00967_43_000000000-A3JHG_1_1113_11294_24024"
#> [143] "M00967_43_000000000-A3JHG_1_1113_13777_28036"
#> [144] "M00967_43_000000000-A3JHG_1_1113_16913_23641"
#> [145] "M00967_43_000000000-A3JHG_1_1113_18037_24127"
#> [146] "M00967_43_000000000-A3JHG_1_1113_18548_4292" 
#> [147] "M00967_43_000000000-A3JHG_1_1113_18644_24962"
#> [148] "M00967_43_000000000-A3JHG_1_1113_18958_10702"
#> [149] "M00967_43_000000000-A3JHG_1_1113_19532_4899" 
#> [150] "M00967_43_000000000-A3JHG_1_1113_20770_14973"
#> [151] "M00967_43_000000000-A3JHG_1_1113_22254_26020"
#> [152] "M00967_43_000000000-A3JHG_1_1113_5336_24219" 
#> [153] "M00967_43_000000000-A3JHG_1_1113_8139_17897" 
#> [154] "M00967_43_000000000-A3JHG_1_1114_11769_3401" 
#> [155] "M00967_43_000000000-A3JHG_1_1114_12471_18325"
#> [156] "M00967_43_000000000-A3JHG_1_1114_13397_2211" 
#> [157] "M00967_43_000000000-A3JHG_1_1114_14280_26421"
#> [158] "M00967_43_000000000-A3JHG_1_1114_14431_2336" 
#> [159] "M00967_43_000000000-A3JHG_1_1114_19875_9916" 
#> [160] "M00967_43_000000000-A3JHG_1_1114_21222_13060"
#> [161] "M00967_43_000000000-A3JHG_1_1114_21515_19829"
#> [162] "M00967_43_000000000-A3JHG_1_1114_22681_21130"
#> [163] "M00967_43_000000000-A3JHG_1_1114_22804_12324"
#> [164] "M00967_43_000000000-A3JHG_1_1114_23055_10307"
#> [165] "M00967_43_000000000-A3JHG_1_1114_25590_23883"
#> [166] "M00967_43_000000000-A3JHG_1_1114_25962_10907"
#> [167] "M00967_43_000000000-A3JHG_1_1114_26968_12650"
#> [168] "M00967_43_000000000-A3JHG_1_1114_3131_16050" 
#> [169] "M00967_43_000000000-A3JHG_1_1114_5902_6841"  
#> [170] "M00967_43_000000000-A3JHG_1_1114_8676_16145" 
#> [171] "M00967_43_000000000-A3JHG_1_2101_13105_15429"
#> [172] "M00967_43_000000000-A3JHG_1_2101_13595_15201"
#> [173] "M00967_43_000000000-A3JHG_1_2101_15190_13450"
#> [174] "M00967_43_000000000-A3JHG_1_2101_16955_14005"
#> [175] "M00967_43_000000000-A3JHG_1_2101_21575_17061"
#> [176] "M00967_43_000000000-A3JHG_1_2101_22400_13416"
#> [177] "M00967_43_000000000-A3JHG_1_2101_22942_18205"
#> [178] "M00967_43_000000000-A3JHG_1_2101_27541_14282"
#> [179] "M00967_43_000000000-A3JHG_1_2101_6377_14009" 
#> [180] "M00967_43_000000000-A3JHG_1_2102_12263_16584"
#> [181] "M00967_43_000000000-A3JHG_1_2102_12696_27241"
#> [182] "M00967_43_000000000-A3JHG_1_2102_15356_3054" 
#> [183] "M00967_43_000000000-A3JHG_1_2102_15542_26928"
#> [184] "M00967_43_000000000-A3JHG_1_2102_15692_18713"
#> [185] "M00967_43_000000000-A3JHG_1_2102_15728_25058"
#> [186] "M00967_43_000000000-A3JHG_1_2102_17019_10818"
#> [187] "M00967_43_000000000-A3JHG_1_2102_27520_15637"
#> [188] "M00967_43_000000000-A3JHG_1_2102_6705_15205" 
#> [189] "M00967_43_000000000-A3JHG_1_2103_10155_24625"
#> [190] "M00967_43_000000000-A3JHG_1_2103_11533_13455"
#> [191] "M00967_43_000000000-A3JHG_1_2103_11635_15432"
#> [192] "M00967_43_000000000-A3JHG_1_2103_12830_18492"
#> [193] "M00967_43_000000000-A3JHG_1_2103_15275_25227"
#> [194] "M00967_43_000000000-A3JHG_1_2103_17221_27256"
#> [195] "M00967_43_000000000-A3JHG_1_2103_20332_24524"
#> [196] "M00967_43_000000000-A3JHG_1_2103_21597_21834"
#> [197] "M00967_43_000000000-A3JHG_1_2103_22029_15973"
#> [198] "M00967_43_000000000-A3JHG_1_2103_22514_6355" 
#> [199] "M00967_43_000000000-A3JHG_1_2103_25452_6018" 
#> [200] "M00967_43_000000000-A3JHG_1_2103_25809_24518"
#> [201] "M00967_43_000000000-A3JHG_1_2104_12407_10383"
#> [202] "M00967_43_000000000-A3JHG_1_2104_12505_3424" 
#> [203] "M00967_43_000000000-A3JHG_1_2104_12560_11750"
#> [204] "M00967_43_000000000-A3JHG_1_2104_13679_3191" 
#> [205] "M00967_43_000000000-A3JHG_1_2104_15263_18242"
#> [206] "M00967_43_000000000-A3JHG_1_2104_17426_10838"
#> [207] "M00967_43_000000000-A3JHG_1_2104_20685_7431" 
#> [208] "M00967_43_000000000-A3JHG_1_2104_22645_18534"
#> [209] "M00967_43_000000000-A3JHG_1_2104_25418_11257"
#> [210] "M00967_43_000000000-A3JHG_1_2104_26311_10309"
#> [211] "M00967_43_000000000-A3JHG_1_2105_10725_6078" 
#> [212] "M00967_43_000000000-A3JHG_1_2105_12878_11098"
#> [213] "M00967_43_000000000-A3JHG_1_2105_17320_11249"
#> [214] "M00967_43_000000000-A3JHG_1_2105_17957_13444"
#> [215] "M00967_43_000000000-A3JHG_1_2105_20695_7710" 
#> [216] "M00967_43_000000000-A3JHG_1_2105_26628_10360"
#> [217] "M00967_43_000000000-A3JHG_1_2105_8155_16465" 
#> [218] "M00967_43_000000000-A3JHG_1_2106_10758_21336"
#> [219] "M00967_43_000000000-A3JHG_1_2106_11005_4817" 
#> [220] "M00967_43_000000000-A3JHG_1_2106_12340_14438"
#> [221] "M00967_43_000000000-A3JHG_1_2106_12968_23126"
#> [222] "M00967_43_000000000-A3JHG_1_2106_14236_13388"
#> [223] "M00967_43_000000000-A3JHG_1_2106_14305_11884"
#> [224] "M00967_43_000000000-A3JHG_1_2106_17516_3721" 
#> [225] "M00967_43_000000000-A3JHG_1_2106_19750_8249" 
#> [226] "M00967_43_000000000-A3JHG_1_2106_20954_22445"
#> [227] "M00967_43_000000000-A3JHG_1_2106_23199_7043" 
#> [228] "M00967_43_000000000-A3JHG_1_2106_25519_17332"
#> [229] "M00967_43_000000000-A3JHG_1_2106_27480_14501"
#> [230] "M00967_43_000000000-A3JHG_1_2106_8736_3855"  
#> [231] "M00967_43_000000000-A3JHG_1_2107_11601_20781"
#> [232] "M00967_43_000000000-A3JHG_1_2107_12056_5129" 
#> [233] "M00967_43_000000000-A3JHG_1_2107_12502_13322"
#> [234] "M00967_43_000000000-A3JHG_1_2107_18723_18535"
#> [235] "M00967_43_000000000-A3JHG_1_2107_25731_13155"
#> [236] "M00967_43_000000000-A3JHG_1_2107_26555_14955"
#> [237] "M00967_43_000000000-A3JHG_1_2107_26817_14749"
#> [238] "M00967_43_000000000-A3JHG_1_2107_6956_16009" 
#> [239] "M00967_43_000000000-A3JHG_1_2107_8240_20277" 
#> [240] "M00967_43_000000000-A3JHG_1_2108_18470_11561"
#> [241] "M00967_43_000000000-A3JHG_1_2108_19309_25091"
#> [242] "M00967_43_000000000-A3JHG_1_2108_19840_26409"
#> [243] "M00967_43_000000000-A3JHG_1_2108_20440_11777"
#> [244] "M00967_43_000000000-A3JHG_1_2108_21552_11279"
#> [245] "M00967_43_000000000-A3JHG_1_2108_22213_11738"
#> [246] "M00967_43_000000000-A3JHG_1_2108_22507_11051"
#> [247] "M00967_43_000000000-A3JHG_1_2108_24812_22312"
#> [248] "M00967_43_000000000-A3JHG_1_2108_25681_15613"
#> [249] "M00967_43_000000000-A3JHG_1_2108_27907_12796"
#> [250] "M00967_43_000000000-A3JHG_1_2108_4168_17142" 
#> [251] "M00967_43_000000000-A3JHG_1_2108_5077_12621" 
#> [252] "M00967_43_000000000-A3JHG_1_2108_6318_9727"  
#> [253] "M00967_43_000000000-A3JHG_1_2108_7674_17990" 
#> [254] "M00967_43_000000000-A3JHG_1_2108_8189_19210" 
#> [255] "M00967_43_000000000-A3JHG_1_2108_8220_23281" 
#> [256] "M00967_43_000000000-A3JHG_1_2109_13206_18588"
#> [257] "M00967_43_000000000-A3JHG_1_2109_17345_6668" 
#> [258] "M00967_43_000000000-A3JHG_1_2109_21549_9949" 
#> [259] "M00967_43_000000000-A3JHG_1_2109_9943_10858" 
#> [260] "M00967_43_000000000-A3JHG_1_2110_10107_20150"
#> [261] "M00967_43_000000000-A3JHG_1_2110_10606_15947"
#> [262] "M00967_43_000000000-A3JHG_1_2110_16815_14755"
#> [263] "M00967_43_000000000-A3JHG_1_2110_18978_25958"
#> [264] "M00967_43_000000000-A3JHG_1_2110_20054_3531" 
#> [265] "M00967_43_000000000-A3JHG_1_2110_21191_16122"
#> [266] "M00967_43_000000000-A3JHG_1_2110_21271_27854"
#> [267] "M00967_43_000000000-A3JHG_1_2110_23478_21034"
#> [268] "M00967_43_000000000-A3JHG_1_2110_23584_15948"
#> [269] "M00967_43_000000000-A3JHG_1_2110_27497_19046"
#> [270] "M00967_43_000000000-A3JHG_1_2110_4035_19768" 
#> [271] "M00967_43_000000000-A3JHG_1_2110_6116_13104" 
#> [272] "M00967_43_000000000-A3JHG_1_2110_8957_20034" 
#> [273] "M00967_43_000000000-A3JHG_1_2110_9069_16610" 
#> [274] "M00967_43_000000000-A3JHG_1_2110_9142_17889" 
#> [275] "M00967_43_000000000-A3JHG_1_2111_10373_11883"
#> [276] "M00967_43_000000000-A3JHG_1_2111_13741_6972" 
#> [277] "M00967_43_000000000-A3JHG_1_2111_14452_16768"
#> [278] "M00967_43_000000000-A3JHG_1_2111_15748_5366" 
#> [279] "M00967_43_000000000-A3JHG_1_2111_16376_7553" 
#> [280] "M00967_43_000000000-A3JHG_1_2111_17223_18856"
#> [281] "M00967_43_000000000-A3JHG_1_2111_22278_7396" 
#> [282] "M00967_43_000000000-A3JHG_1_2111_29227_17392"
#> [283] "M00967_43_000000000-A3JHG_1_2111_5664_6671"  
#> [284] "M00967_43_000000000-A3JHG_1_2111_6782_8678"  
#> [285] "M00967_43_000000000-A3JHG_1_2111_6887_10490" 
#> [286] "M00967_43_000000000-A3JHG_1_2111_8444_24070" 
#> [287] "M00967_43_000000000-A3JHG_1_2112_12596_22917"
#> [288] "M00967_43_000000000-A3JHG_1_2112_15899_17037"
#> [289] "M00967_43_000000000-A3JHG_1_2112_19543_24449"
#> [290] "M00967_43_000000000-A3JHG_1_2112_21752_23067"
#> [291] "M00967_43_000000000-A3JHG_1_2112_23618_22699"
#> [292] "M00967_43_000000000-A3JHG_1_2112_26066_19401"
#> [293] "M00967_43_000000000-A3JHG_1_2112_26773_14017"
#> [294] "M00967_43_000000000-A3JHG_1_2112_6716_23007" 
#> [295] "M00967_43_000000000-A3JHG_1_2112_8202_24000" 
#> [296] "M00967_43_000000000-A3JHG_1_2113_12406_27672"
#> [297] "M00967_43_000000000-A3JHG_1_2113_14487_9130" 
#> [298] "M00967_43_000000000-A3JHG_1_2113_18397_7516" 
#> [299] "M00967_43_000000000-A3JHG_1_2113_20030_3857" 
#> [300] "M00967_43_000000000-A3JHG_1_2113_21028_10433"
#> [301] "M00967_43_000000000-A3JHG_1_2113_21352_15666"
#> [302] "M00967_43_000000000-A3JHG_1_2113_23466_13467"
#> [303] "M00967_43_000000000-A3JHG_1_2113_7021_18070" 
#> [304] "M00967_43_000000000-A3JHG_1_2113_9708_9804"  
#> [305] "M00967_43_000000000-A3JHG_1_2114_10492_14599"
#> [306] "M00967_43_000000000-A3JHG_1_2114_16946_3813" 
#> [307] "M00967_43_000000000-A3JHG_1_2114_27537_17929"
#> [308] "M00967_43_000000000-A3JHG_1_2114_4298_19710" 
#> [309] "M00967_43_000000000-A3JHG_1_2114_6259_8307"  

#' # To get the names of the sequences unique to sample 'F3D0'
miseq$names(type = "sequence", samples = c("F3D0"), distinct = TRUE)
#>   [1] "M00967_43_000000000-A3JHG_1_1101_10331_23332"
#>   [2] "M00967_43_000000000-A3JHG_1_1101_11035_15765"
#>   [3] "M00967_43_000000000-A3JHG_1_1101_14364_8401" 
#>   [4] "M00967_43_000000000-A3JHG_1_1101_18346_24737"
#>   [5] "M00967_43_000000000-A3JHG_1_1101_9620_19745" 
#>   [6] "M00967_43_000000000-A3JHG_1_1102_11115_19075"
#>   [7] "M00967_43_000000000-A3JHG_1_1102_14905_6542" 
#>   [8] "M00967_43_000000000-A3JHG_1_1102_15058_11924"
#>   [9] "M00967_43_000000000-A3JHG_1_1102_16785_23687"
#>  [10] "M00967_43_000000000-A3JHG_1_1102_2114_15227" 
#>  [11] "M00967_43_000000000-A3JHG_1_1103_13966_3813" 
#>  [12] "M00967_43_000000000-A3JHG_1_1103_22490_21890"
#>  [13] "M00967_43_000000000-A3JHG_1_1103_3864_17599" 
#>  [14] "M00967_43_000000000-A3JHG_1_1104_10834_8009" 
#>  [15] "M00967_43_000000000-A3JHG_1_1105_11964_4686" 
#>  [16] "M00967_43_000000000-A3JHG_1_1105_14547_12843"
#>  [17] "M00967_43_000000000-A3JHG_1_1105_9179_20196" 
#>  [18] "M00967_43_000000000-A3JHG_1_1106_11240_22282"
#>  [19] "M00967_43_000000000-A3JHG_1_1106_12231_13452"
#>  [20] "M00967_43_000000000-A3JHG_1_1106_26087_11342"
#>  [21] "M00967_43_000000000-A3JHG_1_1106_26888_10548"
#>  [22] "M00967_43_000000000-A3JHG_1_1106_9558_4409"  
#>  [23] "M00967_43_000000000-A3JHG_1_1107_18662_2044" 
#>  [24] "M00967_43_000000000-A3JHG_1_1107_21038_8569" 
#>  [25] "M00967_43_000000000-A3JHG_1_1108_12037_3053" 
#>  [26] "M00967_43_000000000-A3JHG_1_1108_18996_9436" 
#>  [27] "M00967_43_000000000-A3JHG_1_1108_19065_7362" 
#>  [28] "M00967_43_000000000-A3JHG_1_1109_14359_12484"
#>  [29] "M00967_43_000000000-A3JHG_1_1109_24274_5733" 
#>  [30] "M00967_43_000000000-A3JHG_1_1109_7971_24259" 
#>  [31] "M00967_43_000000000-A3JHG_1_1109_8014_5010"  
#>  [32] "M00967_43_000000000-A3JHG_1_1110_13719_21639"
#>  [33] "M00967_43_000000000-A3JHG_1_1110_14443_28339"
#>  [34] "M00967_43_000000000-A3JHG_1_1111_17002_7999" 
#>  [35] "M00967_43_000000000-A3JHG_1_1111_8373_10700" 
#>  [36] "M00967_43_000000000-A3JHG_1_1112_10194_20788"
#>  [37] "M00967_43_000000000-A3JHG_1_1112_11343_21537"
#>  [38] "M00967_43_000000000-A3JHG_1_1112_12844_3769" 
#>  [39] "M00967_43_000000000-A3JHG_1_1112_20884_13058"
#>  [40] "M00967_43_000000000-A3JHG_1_1113_16913_23641"
#>  [41] "M00967_43_000000000-A3JHG_1_1113_18548_4292" 
#>  [42] "M00967_43_000000000-A3JHG_1_1113_18958_10702"
#>  [43] "M00967_43_000000000-A3JHG_1_1114_12471_18325"
#>  [44] "M00967_43_000000000-A3JHG_1_1114_13397_2211" 
#>  [45] "M00967_43_000000000-A3JHG_1_1114_14431_2336" 
#>  [46] "M00967_43_000000000-A3JHG_1_1114_21222_13060"
#>  [47] "M00967_43_000000000-A3JHG_1_1114_22804_12324"
#>  [48] "M00967_43_000000000-A3JHG_1_1114_23055_10307"
#>  [49] "M00967_43_000000000-A3JHG_1_1114_26968_12650"
#>  [50] "M00967_43_000000000-A3JHG_1_1114_3131_16050" 
#>  [51] "M00967_43_000000000-A3JHG_1_1114_5902_6841"  
#>  [52] "M00967_43_000000000-A3JHG_1_1114_8676_16145" 
#>  [53] "M00967_43_000000000-A3JHG_1_2101_16955_14005"
#>  [54] "M00967_43_000000000-A3JHG_1_2101_22942_18205"
#>  [55] "M00967_43_000000000-A3JHG_1_2102_12263_16584"
#>  [56] "M00967_43_000000000-A3JHG_1_2102_12696_27241"
#>  [57] "M00967_43_000000000-A3JHG_1_2102_17019_10818"
#>  [58] "M00967_43_000000000-A3JHG_1_2102_6705_15205" 
#>  [59] "M00967_43_000000000-A3JHG_1_2103_17221_27256"
#>  [60] "M00967_43_000000000-A3JHG_1_2103_22029_15973"
#>  [61] "M00967_43_000000000-A3JHG_1_2103_22514_6355" 
#>  [62] "M00967_43_000000000-A3JHG_1_2103_25452_6018" 
#>  [63] "M00967_43_000000000-A3JHG_1_2104_12505_3424" 
#>  [64] "M00967_43_000000000-A3JHG_1_2105_12878_11098"
#>  [65] "M00967_43_000000000-A3JHG_1_2105_17320_11249"
#>  [66] "M00967_43_000000000-A3JHG_1_2105_26628_10360"
#>  [67] "M00967_43_000000000-A3JHG_1_2105_8155_16465" 
#>  [68] "M00967_43_000000000-A3JHG_1_2106_10758_21336"
#>  [69] "M00967_43_000000000-A3JHG_1_2106_14305_11884"
#>  [70] "M00967_43_000000000-A3JHG_1_2106_19750_8249" 
#>  [71] "M00967_43_000000000-A3JHG_1_2106_27480_14501"
#>  [72] "M00967_43_000000000-A3JHG_1_2107_25731_13155"
#>  [73] "M00967_43_000000000-A3JHG_1_2107_26555_14955"
#>  [74] "M00967_43_000000000-A3JHG_1_2107_26817_14749"
#>  [75] "M00967_43_000000000-A3JHG_1_2107_8240_20277" 
#>  [76] "M00967_43_000000000-A3JHG_1_2108_19840_26409"
#>  [77] "M00967_43_000000000-A3JHG_1_2108_5077_12621" 
#>  [78] "M00967_43_000000000-A3JHG_1_2108_6318_9727"  
#>  [79] "M00967_43_000000000-A3JHG_1_2108_8220_23281" 
#>  [80] "M00967_43_000000000-A3JHG_1_2109_17345_6668" 
#>  [81] "M00967_43_000000000-A3JHG_1_2109_21549_9949" 
#>  [82] "M00967_43_000000000-A3JHG_1_2110_10107_20150"
#>  [83] "M00967_43_000000000-A3JHG_1_2110_10606_15947"
#>  [84] "M00967_43_000000000-A3JHG_1_2110_16815_14755"
#>  [85] "M00967_43_000000000-A3JHG_1_2110_18978_25958"
#>  [86] "M00967_43_000000000-A3JHG_1_2110_20054_3531" 
#>  [87] "M00967_43_000000000-A3JHG_1_2110_23478_21034"
#>  [88] "M00967_43_000000000-A3JHG_1_2110_9069_16610" 
#>  [89] "M00967_43_000000000-A3JHG_1_2111_17223_18856"
#>  [90] "M00967_43_000000000-A3JHG_1_2111_22278_7396" 
#>  [91] "M00967_43_000000000-A3JHG_1_2111_29227_17392"
#>  [92] "M00967_43_000000000-A3JHG_1_2111_6782_8678"  
#>  [93] "M00967_43_000000000-A3JHG_1_2111_6887_10490" 
#>  [94] "M00967_43_000000000-A3JHG_1_2111_8444_24070" 
#>  [95] "M00967_43_000000000-A3JHG_1_2112_15899_17037"
#>  [96] "M00967_43_000000000-A3JHG_1_2112_26066_19401"
#>  [97] "M00967_43_000000000-A3JHG_1_2113_12406_27672"
#>  [98] "M00967_43_000000000-A3JHG_1_2113_18397_7516" 
#>  [99] "M00967_43_000000000-A3JHG_1_2113_9708_9804"  
#> [100] "M00967_43_000000000-A3JHG_1_2114_27537_17929"
#> [101] "M00967_43_000000000-A3JHG_1_2114_4298_19710" 

# To get the names of the samples
miseq$names(type = "sample")
#>  [1] "F3D0"   "F3D1"   "F3D141" "F3D142" "F3D143" "F3D144" "F3D145" "F3D146"
#>  [9] "F3D147" "F3D148" "F3D149" "F3D150" "F3D2"   "F3D3"   "F3D5"   "F3D6"  
#> [17] "F3D7"   "F3D8"   "F3D9"  

# To get the names of the treatments
miseq$names(type = "treatment")
#> [1] "Early" "Late" 

# To get the names of the bins
miseq$names(type = "bin")
#>   [1] "Otu001" "Otu002" "Otu003" "Otu004" "Otu005" "Otu006" "Otu007" "Otu008"
#>   [9] "Otu009" "Otu010" "Otu011" "Otu012" "Otu013" "Otu014" "Otu015" "Otu016"
#>  [17] "Otu017" "Otu018" "Otu019" "Otu020" "Otu021" "Otu022" "Otu023" "Otu024"
#>  [25] "Otu025" "Otu026" "Otu027" "Otu028" "Otu029" "Otu030" "Otu031" "Otu032"
#>  [33] "Otu033" "Otu034" "Otu035" "Otu036" "Otu037" "Otu038" "Otu039" "Otu040"
#>  [41] "Otu041" "Otu042" "Otu043" "Otu044" "Otu045" "Otu046" "Otu047" "Otu048"
#>  [49] "Otu049" "Otu050" "Otu051" "Otu052" "Otu053" "Otu054" "Otu055" "Otu056"
#>  [57] "Otu057" "Otu058" "Otu059" "Otu060" "Otu061" "Otu062" "Otu063" "Otu064"
#>  [65] "Otu065" "Otu066" "Otu067" "Otu068" "Otu069" "Otu070" "Otu071" "Otu072"
#>  [73] "Otu073" "Otu074" "Otu075" "Otu076" "Otu077" "Otu078" "Otu079" "Otu080"
#>  [81] "Otu081" "Otu082" "Otu083" "Otu084" "Otu085" "Otu086" "Otu087" "Otu088"
#>  [89] "Otu089" "Otu090" "Otu091" "Otu092" "Otu093" "Otu094" "Otu095" "Otu096"
#>  [97] "Otu097" "Otu098" "Otu099" "Otu100" "Otu101" "Otu102" "Otu103" "Otu104"
#> [105] "Otu105" "Otu106" "Otu107" "Otu108" "Otu109" "Otu110" "Otu111" "Otu112"
#> [113] "Otu113" "Otu114" "Otu115" "Otu116" "Otu117" "Otu118" "Otu119" "Otu120"
#> [121] "Otu121" "Otu122" "Otu123" "Otu124" "Otu125" "Otu126" "Otu127" "Otu128"
#> [129] "Otu129" "Otu130" "Otu131" "Otu132" "Otu133" "Otu134" "Otu135" "Otu136"
#> [137] "Otu137" "Otu138" "Otu139" "Otu140" "Otu141" "Otu142" "Otu143" "Otu144"
#> [145] "Otu145" "Otu146" "Otu147" "Otu148" "Otu149" "Otu150" "Otu151" "Otu152"
#> [153] "Otu153" "Otu154" "Otu155" "Otu156" "Otu157" "Otu158" "Otu159" "Otu160"
#> [161] "Otu161" "Otu162" "Otu163" "Otu164" "Otu165" "Otu166" "Otu167" "Otu168"
#> [169] "Otu169" "Otu170" "Otu171" "Otu172" "Otu173" "Otu174" "Otu175" "Otu176"
#> [177] "Otu177" "Otu178" "Otu179" "Otu180" "Otu181" "Otu182" "Otu183" "Otu184"
#> [185] "Otu185" "Otu186" "Otu187" "Otu188" "Otu189" "Otu190" "Otu191" "Otu192"
#> [193] "Otu193" "Otu194" "Otu195" "Otu196" "Otu197" "Otu198" "Otu199" "Otu200"
#> [201] "Otu201" "Otu202" "Otu203" "Otu204" "Otu205" "Otu206" "Otu207" "Otu208"
#> [209] "Otu209" "Otu210" "Otu211" "Otu212" "Otu213" "Otu214" "Otu215" "Otu216"
#> [217] "Otu217" "Otu218" "Otu219" "Otu220" "Otu221" "Otu222" "Otu223" "Otu224"
#> [225] "Otu225" "Otu226" "Otu227" "Otu228" "Otu229" "Otu230" "Otu231" "Otu232"
#> [233] "Otu233" "Otu234" "Otu235" "Otu236" "Otu237" "Otu238" "Otu239" "Otu240"
#> [241] "Otu241" "Otu242" "Otu243" "Otu244" "Otu245" "Otu246" "Otu247" "Otu248"
#> [249] "Otu249" "Otu250" "Otu251" "Otu252" "Otu253" "Otu254" "Otu255" "Otu256"
#> [257] "Otu257" "Otu258" "Otu259" "Otu260" "Otu261" "Otu262" "Otu263" "Otu264"
#> [265] "Otu265" "Otu266" "Otu267" "Otu268" "Otu269" "Otu270" "Otu271" "Otu272"
#> [273] "Otu273" "Otu274" "Otu275" "Otu276" "Otu277" "Otu278" "Otu279" "Otu280"
#> [281] "Otu281" "Otu282" "Otu283" "Otu284" "Otu285" "Otu286" "Otu287" "Otu288"
#> [289] "Otu289" "Otu290" "Otu291" "Otu292" "Otu293" "Otu294" "Otu295" "Otu296"
#> [297] "Otu297" "Otu298" "Otu299" "Otu300" "Otu301" "Otu302" "Otu303" "Otu304"
#> [305] "Otu305" "Otu306" "Otu307" "Otu308" "Otu309" "Otu310" "Otu311" "Otu312"
#> [313] "Otu313" "Otu314" "Otu315" "Otu316" "Otu317" "Otu318" "Otu319" "Otu320"
#> [321] "Otu321" "Otu322" "Otu323" "Otu324" "Otu325" "Otu326" "Otu327" "Otu328"
#> [329] "Otu329" "Otu330" "Otu331" "Otu332" "Otu333" "Otu334" "Otu335" "Otu336"
#> [337] "Otu337" "Otu338" "Otu339" "Otu340" "Otu341" "Otu342" "Otu343" "Otu344"
#> [345] "Otu345" "Otu346" "Otu347" "Otu348" "Otu349" "Otu350" "Otu351" "Otu352"
#> [353] "Otu353" "Otu354" "Otu355" "Otu356" "Otu357" "Otu358" "Otu359" "Otu360"
#> [361] "Otu361" "Otu362" "Otu363" "Otu364" "Otu365" "Otu366" "Otu367" "Otu368"
#> [369] "Otu369" "Otu370" "Otu371" "Otu372" "Otu373" "Otu374" "Otu375" "Otu376"
#> [377] "Otu377" "Otu378" "Otu379" "Otu380" "Otu381" "Otu382" "Otu383" "Otu384"
#> [385] "Otu385" "Otu386" "Otu387" "Otu388" "Otu389" "Otu390" "Otu391" "Otu392"
#> [393] "Otu393" "Otu394" "Otu395" "Otu396" "Otu397" "Otu398" "Otu399" "Otu400"
#> [401] "Otu401" "Otu402" "Otu403" "Otu404" "Otu405" "Otu406" "Otu407" "Otu408"
#> [409] "Otu409" "Otu410" "Otu411" "Otu412" "Otu413" "Otu414" "Otu415" "Otu416"
#> [417] "Otu417" "Otu418" "Otu419" "Otu420" "Otu421" "Otu422" "Otu423" "Otu424"
#> [425] "Otu425" "Otu426" "Otu427" "Otu428" "Otu429" "Otu430" "Otu431" "Otu432"
#> [433] "Otu433" "Otu434" "Otu435" "Otu436" "Otu437" "Otu438" "Otu439" "Otu440"
#> [441] "Otu441" "Otu442" "Otu443" "Otu444" "Otu445" "Otu446" "Otu447" "Otu448"
#> [449] "Otu449" "Otu450" "Otu451" "Otu452" "Otu453" "Otu454" "Otu455" "Otu456"
#> [457] "Otu457" "Otu458" "Otu459" "Otu460" "Otu461" "Otu462" "Otu463" "Otu464"
#> [465] "Otu465" "Otu466" "Otu467" "Otu468" "Otu469" "Otu470" "Otu471" "Otu472"
#> [473] "Otu473" "Otu474" "Otu475" "Otu476" "Otu477" "Otu478" "Otu479" "Otu480"
#> [481] "Otu481" "Otu482" "Otu483" "Otu484" "Otu485" "Otu486" "Otu487" "Otu488"
#> [489] "Otu489" "Otu490" "Otu491" "Otu492" "Otu493" "Otu494" "Otu495" "Otu496"
#> [497] "Otu497" "Otu498" "Otu499" "Otu500" "Otu501" "Otu502" "Otu503" "Otu504"
#> [505] "Otu505" "Otu506" "Otu507" "Otu508" "Otu509" "Otu510" "Otu511" "Otu512"
#> [513] "Otu513" "Otu514" "Otu515" "Otu516" "Otu517" "Otu518" "Otu519" "Otu520"
#> [521] "Otu521" "Otu522" "Otu523" "Otu524" "Otu525" "Otu526" "Otu527" "Otu528"
#> [529] "Otu529" "Otu530" "Otu531"

# To get the names of the bins that are unique to 'F3D0'
miseq$names(type = "bin", samples = c("F3D0"), distinct = TRUE)
#>  [1] "Otu330" "Otu339" "Otu341" "Otu345" "Otu347" "Otu354" "Otu364" "Otu431"
#>  [9] "Otu466" "Otu469" "Otu470" "Otu491" "Otu493" "Otu529"

# To get the names of the bins that include sequences from 'F3D0'
miseq$names(type = "bin", samples = c("F3D0"), distinct = FALSE)
#>   [1] "Otu001" "Otu002" "Otu003" "Otu004" "Otu005" "Otu006" "Otu007" "Otu008"
#>   [9] "Otu009" "Otu010" "Otu011" "Otu012" "Otu013" "Otu014" "Otu015" "Otu016"
#>  [17] "Otu017" "Otu018" "Otu019" "Otu020" "Otu021" "Otu022" "Otu023" "Otu024"
#>  [25] "Otu025" "Otu026" "Otu027" "Otu028" "Otu029" "Otu030" "Otu031" "Otu032"
#>  [33] "Otu033" "Otu034" "Otu035" "Otu036" "Otu037" "Otu038" "Otu039" "Otu040"
#>  [41] "Otu041" "Otu042" "Otu043" "Otu044" "Otu045" "Otu046" "Otu047" "Otu048"
#>  [49] "Otu049" "Otu050" "Otu051" "Otu052" "Otu053" "Otu054" "Otu056" "Otu057"
#>  [57] "Otu060" "Otu061" "Otu062" "Otu063" "Otu064" "Otu065" "Otu066" "Otu067"
#>  [65] "Otu068" "Otu069" "Otu070" "Otu071" "Otu072" "Otu073" "Otu074" "Otu075"
#>  [73] "Otu076" "Otu078" "Otu079" "Otu082" "Otu083" "Otu084" "Otu085" "Otu087"
#>  [81] "Otu089" "Otu090" "Otu091" "Otu092" "Otu093" "Otu094" "Otu095" "Otu096"
#>  [89] "Otu097" "Otu098" "Otu101" "Otu102" "Otu103" "Otu104" "Otu105" "Otu106"
#>  [97] "Otu107" "Otu108" "Otu109" "Otu111" "Otu112" "Otu113" "Otu114" "Otu115"
#> [105] "Otu117" "Otu118" "Otu119" "Otu120" "Otu122" "Otu123" "Otu124" "Otu125"
#> [113] "Otu127" "Otu129" "Otu130" "Otu132" "Otu136" "Otu137" "Otu138" "Otu141"
#> [121] "Otu143" "Otu144" "Otu145" "Otu146" "Otu147" "Otu148" "Otu150" "Otu151"
#> [129] "Otu154" "Otu156" "Otu157" "Otu160" "Otu161" "Otu162" "Otu163" "Otu164"
#> [137] "Otu166" "Otu168" "Otu169" "Otu170" "Otu172" "Otu174" "Otu175" "Otu176"
#> [145] "Otu178" "Otu182" "Otu183" "Otu186" "Otu187" "Otu188" "Otu190" "Otu196"
#> [153] "Otu197" "Otu202" "Otu203" "Otu206" "Otu207" "Otu208" "Otu211" "Otu215"
#> [161] "Otu216" "Otu220" "Otu222" "Otu224" "Otu225" "Otu227" "Otu237" "Otu243"
#> [169] "Otu256" "Otu259" "Otu275" "Otu277" "Otu287" "Otu290" "Otu300" "Otu310"
#> [177] "Otu312" "Otu330" "Otu339" "Otu341" "Otu345" "Otu347" "Otu354" "Otu364"
#> [185] "Otu431" "Otu466" "Otu469" "Otu470" "Otu491" "Otu493" "Otu529"

# To get the names of the reports
miseq$names(type = "report")
#> [1] "contigs_report" "metadata"      


## ------------------------------------------------
## Method `strollur$report()`
## ------------------------------------------------


miseq <- load_dataset(strollur_example("miseq_sop.rds"))

# To get the FASTA data

miseq$report(type = "fasta") |> head(n = 5)
#>                                  sequence_name
#> 1  M00967_43_000000000-A3JHG_1_1101_10133_8460
#> 2 M00967_43_000000000-A3JHG_1_1101_10331_23332
#> 3 M00967_43_000000000-A3JHG_1_1101_10382_22128
#> 4 M00967_43_000000000-A3JHG_1_1101_11035_15765
#> 5 M00967_43_000000000-A3JHG_1_1101_11348_22601
#>                                                                                                                                                                                                                                                                                                                                                                                  sequence
#> 1 TAC--GT-AG-GGG--GCA-A-G-C-G-T-T--AT-C-CGG-AT--TT-A-C-T--GG-GT--GT-A-AA-GG-GA-GC-G-TA-GGC-G-G-C-CA-T-G-C-AA-G-T-C-A-G-A-A-G--TG-A-AA-AC-C-C-GG-GG--CT-C-AA-C---C-C-TGG-G-AGT-G-C-TTTT-GAAAC-TG-T-GCGGC-TAGA-GT-GT-CG-GA-G-G---GG-T-A-AGTGGAATTCCTAGTGT-AGCGGT-GAAATGCGTAG-AT-A-TT-AG-GA-GG-AACACCAGT-GGCGAAGGCG------GCTTA-CTG-G-AC-GA-T-G-ACTGACG-CTGA-GGCT-CGAAA-GCG-TGGGG-AGC-AAACAGG
#> 2 TAC--GG-AG-GAT--GCG-A-G-C-G-T-T--AT-C-CGG-AT--TT-A-T-T--GG-GT--TT-A-AA-GG-GA-GC-G-CA-GGC-G-G-C-AT-G-G-C-AA-G-T-C-A-G-A-T-G--TG-A-AA-GC-C-C-GG-GG--CT-C-AA-C-C-C-C-G-G-G-ACT-G-C-ATTT-GAAAC-TG-C-CAGGC-TAGA-GT-GT-CG-GA-G-A---GG-C-A-AGTGGAATTCCTAGTGT-AGCGGT-GAAATGCGTAG-AT-A-TT-AG-GA-GG-AACACCAGT-GGCGGAGGCG------GCTTA-CTG-G-AC-GG-T-C-ACTGACG-CTGA-GGCT-CGAAA-GCG-TGGGG-AGC-AAACAGG
#> 3 TAC--GT-AG-GTA--GCA-A-G-C-G-T-T--GT-C-CGG-AT--TT-A-C-T--GG-GT--GT-A-AA-GG-GC-GT-G-TA-GCC-G-G-G-CT-T-A-C-AA-G-T-C-A-G-A-T-G--TG-A-AA-TC-C-G-GG-GG--CT-C-AA-C-C-C-C-C-G-A-ACT-G-C-ATTT-GAAAC-TG-T-AGGTC-TTGA-GT-AT-CG-GA-G-A---GG-C-A-GGCGGAATTCCTAGTGT-AGCGGT-GAAATGCGTAG-AT-A-TT-AG-GA-GG-AACACCAGT-GGCGAAGGCG------GCCTG-CTG-G-AC-GA-C-A-ACTGACG-GTGA-GGCG-CGAAA-GCG-TGGGG-AGC-AAACAGG
#> 4 TAC--GT-AG-GGG--GCA-A-G-C-G-T-T--AT-C-CGG-AT--TT-A-C-T--GG-GT--GT-A-AA-GG-GG-GC-G-TA-GAC-G-G-C-AG-T-G-C-AA-G-T-C-T-G-A-A-G--TG-A-AA-GG-C-G-GG-GG--CC-C-AA-C-C-C-C-C-G-A-ACT-G-C-TTTG-GAAAC-TG-T-GCTGC-TGGA-GT-GC-AG-GA-G-A---GG-C-A-AGTGGAATTCCTAGTGT-AGCGGT-GAAATGCGTAG-AT-A-TT-AG-GA-GG-AACACCAG--GGCGAAGGCG------GCTTG-CTG-G-AC-TG-T-A-ACTGACG-TTGA-GGCT-CGAAA-GCG-TGGGG-AGC-AAACAGG
#> 5 TAC--GG-AG-GAT--GCG-A-G-C-G-T-T--AT-C-CGG-AT--TC-A-C-T--GG-GT--GT-A-AA-GG-GA-GC-G-TA-GAC-G-G-C-CA-T-G-C-AA-G-C-C-A-G-G-G-G--TG-A-AA-GC-C-C-GG-GG--CC-C-AA-C-C-C-C-G-G-G-ACT-G-C-CCTT-GGAAC-TG-C-ATGGC-TGGA-GT-GC-GG-GA-G-G---GG-C-A-GGCGGAATTCCTGGTGT-AGCGGT-GAAATGCGTAG-AT-A-TC-AG-GA-GG-AACACCGGC-GGCGAAGGCG------GCCTG-CTG-G-AC-CG-C-G-ACTGACG-TTGA-GGCT-CGAAA-GCG-TGGGG-AGC-GAACAGG

# To get a report about the FASTA data

miseq$report(type = "sequence") |> head(n = 5)
#>                                  sequence_name start end length ambig
#> 1  M00967_43_000000000-A3JHG_1_1101_10133_8460     1 375    253     0
#> 2 M00967_43_000000000-A3JHG_1_1101_10331_23332     1 375    253     0
#> 3 M00967_43_000000000-A3JHG_1_1101_10382_22128     1 375    253     0
#> 4 M00967_43_000000000-A3JHG_1_1101_11035_15765     1 375    252     0
#> 5 M00967_43_000000000-A3JHG_1_1101_11348_22601     1 375    253     0
#>   longest_homopolymer num_n
#> 1                   5     0
#> 2                   4     0
#> 3                   5     0
#> 4                   5     0
#> 5                   4     0

# To get the sequence bin assignments

miseq$report(type = "sequence_bin_assignment", bin_type = "otu") |>
head(n = 5)
#>   otu_id                                       seq_id
#> 1 Otu001 M00967_43_000000000-A3JHG_1_1101_12302_23776
#> 2 Otu001 M00967_43_000000000-A3JHG_1_1101_20262_22075
#> 3 Otu001 M00967_43_000000000-A3JHG_1_1102_18640_14309
#> 4 Otu001  M00967_43_000000000-A3JHG_1_1102_20738_4913
#> 5 Otu001   M00967_43_000000000-A3JHG_1_1102_6774_6343

# To get the sample treatment assignments

miseq$report(type = "sample_assignment") |> head(n = 5)
#>   sample treatment
#> 1   F3D0     Early
#> 2   F3D1     Early
#> 3 F3D141      Late
#> 4 F3D142      Late
#> 5 F3D143      Late

# To get a report about sequence classifications

miseq$report(type = "sequence_taxonomy") |> head(n = 5)
#>                                 sequence_name level        taxonomy confidence
#> 1 M00967_43_000000000-A3JHG_1_1101_10133_8460     1        Bacteria        100
#> 2 M00967_43_000000000-A3JHG_1_1101_10133_8460     2      Firmicutes        100
#> 3 M00967_43_000000000-A3JHG_1_1101_10133_8460     3      Clostridia        100
#> 4 M00967_43_000000000-A3JHG_1_1101_10133_8460     4   Clostridiales        100
#> 5 M00967_43_000000000-A3JHG_1_1101_10133_8460     5 Lachnospiraceae         85

# To get a report about bin classifications for 'otu' data

miseq$report(type = "bin_taxonomy", bin_type = "otu") |> head(n = 5)
#>   bin_name level             taxonomy confidence
#> 1   Otu001     1             Bacteria        100
#> 2   Otu001     2      "Bacteroidetes"        100
#> 3   Otu001     3        "Bacteroidia"        100
#> 4   Otu001     4      "Bacteroidales"        100
#> 5   Otu001     5 "Porphyromonadaceae"        100

# To get the 'otu' bin representative sequences

miseq$report(type = "bin_representative", bin_type = "otu") |>
head(n = 5)
#>   otu_names                          representative_name
#> 1    Otu001 M00967_43_000000000-A3JHG_1_1108_14299_17220
#> 2    Otu002  M00967_43_000000000-A3JHG_1_1106_22705_6123
#> 3    Otu003  M00967_43_000000000-A3JHG_1_1101_15533_5293
#> 4    Otu004 M00967_43_000000000-A3JHG_1_1105_25642_17588
#> 5    Otu005  M00967_43_000000000-A3JHG_1_2102_7041_13746
#>                                                                                                                                                                                                                                                                                                                                                                   representative_sequence
#> 1 TAC--GT-AG-GGG--GCA-A-G-C-G-T-T--AT-C-CGG-AT--TT-A-C-T--GG-GT--GT-A-AA-GG-GA-GC-G-TA-GAC-G-G-C-TG-T-G-C-AA-G-T-C-T-G-A-A-G--TG-A-AA-TG-C-C-GG-GG--CT-C-AA-C-C-C-C-G-G-A-ACT-G-C-TTTG-GAAAC-TG-T-ACAGC-TAGA-GT-GC-AG-GA-G-G---GG-T-G-AGCGGAATTCCTAGTGT-AGCGGT-GAAATGCGTAG-AT-A-TT-AG-GA-GG-AACACCGGT-GGCGAAGGCG------GCTCA-CTG-G-AC-TG-T-A-ACTGACG-TTGA-GGCT-CGAAA-GCG-TGGGG-AGC-AAACAGG
#> 2 TAC--GT-AG-GGG--GCA-A-G-C-G-T-T--AT-C-CGG-AT--TT-A-C-T--GG-GT--GT-A-AA-GG-GA-GC-G-CA-GAC-G-G-C-TG-T-G-C-AA-G-T-C-T-G-G-A-G--TG-A-AA-GG-C-G-GG-GG--CC-C-AA-C-C-C-C-C-G-G-ACT-G-C-TCTG-GAAAC-TG-T-AAAGC-TGGA-GT-GC-AG-GA-G-A---GG-T-A-AGCGGAATTCCTAGTGT-AGCGGT-GAAATGCGTAG-AT-A-TT-AG-GA-GG-AACACCAGT-GGCGAAGGCG------GCTTA-CTG-G-AC-TG-C-A-ACTGACG-TTGA-GGCT-CGAAA-GCG-TGGGT-ATC-GAACAGG
#> 3 TAC--GG-AG-GAT--GCG-A-G-C-G-T-T--AT-C-CGG-AT--TT-A-C-T--GG-GT--GT-A-AA-GG-GA-GC-G-TA-GAC-G-G-C-GA-T-G-C-AA-G-T-C-T-G-A-A-G--TG-A-AA-GG-C-G-GG-GG--CC-C-AA-C-C-C-C-C-G-G-ACT-G-C-TTTG-GAAAC-TG-T-ATAGC-TGGA-GT-GC-AG-GA-G-A---GG-T-A-AGTGGAATTCCTAGTGT-AGCGGT-GAAATGCGTAG-AT-A-TT-AG-GA-GG-AACACCAGT-GGCGAAGGCG------GCTTA-CTG-G-AC-TG-T-A-ACTGACG-TTGA-GGCT-CGAAA-GCG-TGGGG-AGC-AAACAGG
#> 4 TAC--GT-AG-GTG--GCA-A-G-C-G-T-T--AT-C-CGG-AT--TT-A-C-T--GG-GT--GT-A-AA-GG-GC-GT-G-TA-GGC-G-G-G-AC-T-G-C-AA-G-T-C-A-G-A-T-G--TG-A-AA-CC-C-A-TG-GG--CT-C-AA-C-C-C-A-T-G-G-CCT-G-C-ATTT-GAAAC-TG-T-AGTTC-TTGA-GT-GA-TG-GA-G-A---GG-C-A-GGCGGAATTCCGTGTGT-AGCGGT-GAAATGCGTAG-AT-A-TA-CG-GA-GG-AACACCAGT-GGCGAAGGCG------GCCTG-CTG-G-AC-AT-T-A-ACTGACG-CTGA-GGCG-CGAAA-GCG-TGGGG-AGC-AAACAGG
#> 5 TAC--GT-AG-GGG--GCG-A-G-C-G-T-T--AT-C-CGG-AT--TC-A-T-T--GG-GC--GT-A-AA-GC-GC-GC-G-CA-GGC-G-G-A-CT-C-A-T-AA-G-C-G-G-A-G-C-C--TT-T-AA-TC-T-T-GG-GG--CT-T-AA-C-C-T-C-A-A-G-T-C-G-G-GCCC-CGAAC-TG-T-GAGTC-TCGA-GT-GT-GG-TA-G-G---GG-A-A-GGCGGAATTCCCGGTGT-AGCGGT-GGAATGCGCAG-AT-A-TC-GG-GA-AG-AACACCGAT-GGCGAAGGCA------GCCTT-CTG-G-GC-CA-T-C-ACTGACG-CTGA-GGCG-CGAAA-GCT-AGGGG-AGC-AAACAGG

# To get a report about the sequences removed during your analysis:

miseq$report(type = "sequence_scrap")
#> data frame with 0 columns and 0 rows

# To get a report about the "otu" bins removed during your analysis:

miseq$report(type = "bin_scrap", bin_type = "otu")
#> data frame with 0 columns and 0 rows

# To get the metadata associated with your data:

metadata <- miseq$report(type = "metadata") |> head(n = 5)

# To get the resource references associated with your data:

references <- miseq$report(type = "resource_reference")

# To get our custom report containing the contigs assembly data:

miseq$report(type = "contigs_report") |> head(n = 5)
#>                                           Name Length Overlap_Length
#> 1  M00967_43_000000000-A3JHG_1_1101_10133_8460    253            249
#> 2 M00967_43_000000000-A3JHG_1_1101_10331_23332    253            249
#> 3 M00967_43_000000000-A3JHG_1_1101_10382_22128    253            249
#> 4 M00967_43_000000000-A3JHG_1_1101_11035_15765    252            250
#> 5 M00967_43_000000000-A3JHG_1_1101_11348_22601    253            249
#>   Overlap_Start Overlap_End MisMatches Num_Ns Expected_Errors
#> 1             2         251          0      0      0.00207161
#> 2             2         251          0      0      0.00230661
#> 3             2         251          2      0      0.01134540
#> 4             1         251         19      0      0.08720880
#> 5             2         251          2      0      0.00674431


## ------------------------------------------------
## Method `strollur$summary()`
## ------------------------------------------------


miseq <- load_dataset(strollur_example("miseq_sop.rds"))

# To get the summary of your FASTA data
miseq$summary(type = "sequence")
#>             starts ends   nbases ambigs polymers numns    numseqs
#> Minimum:         1  375 249.0000      0 3.000000     0      1.000
#> 2.5%-tile:       1  375 252.0000      0 4.000000     0   2849.075
#> 25%-tile:        1  375 252.0000      0 4.000000     0  28490.750
#> Median:          1  375 253.0000      0 4.000000     0  56981.500
#> 75%-tile:        1  375 253.0000      0 5.000000     0  85472.250
#> 97.5%-tile:      1  375 254.0000      0 6.000000     0 111113.925
#> Maximum:         1  375 256.0000      0 6.000000     0 113963.000
#> Mean:            1  375 252.7406      0 4.496082     0  56981.643
#>             starts ends   nbases ambigs polymers numns    numseqs
#> Minimum:         1  375 249.0000      0 3.000000     0      1.000
#> 2.5%-tile:       1  375 252.0000      0 4.000000     0   2849.075
#> 25%-tile:        1  375 252.0000      0 4.000000     0  28490.750
#> Median:          1  375 253.0000      0 4.000000     0  56981.500
#> 75%-tile:        1  375 253.0000      0 5.000000     0  85472.250
#> 97.5%-tile:      1  375 254.0000      0 6.000000     0 111113.925
#> Maximum:         1  375 256.0000      0 6.000000     0 113963.000
#> Mean:            1  375 252.7406      0 4.496082     0  56981.643
#>             starts ends   nbases ambigs polymers numns    numseqs
#> Minimum:         1  375 249.0000      0 3.000000     0      1.000
#> 2.5%-tile:       1  375 252.0000      0 4.000000     0   2849.075
#> 25%-tile:        1  375 252.0000      0 4.000000     0  28490.750
#> Median:          1  375 253.0000      0 4.000000     0  56981.500
#> 75%-tile:        1  375 253.0000      0 5.000000     0  85472.250
#> 97.5%-tile:      1  375 254.0000      0 6.000000     0 111113.925
#> Maximum:         1  375 256.0000      0 6.000000     0 113963.000
#> Mean:            1  375 252.7406      0 4.496082     0  56981.643

# summarize contigs_report
miseq$summary(type = "report", report_type = "contigs_report")
#>               Length Overlap_Length Overlap_Start Overlap_End MisMatches Num_Ns
#> Minimum:    250.0000       232.0000      0.000000    248.0000   0.000000      0
#> 2.5%-tile:  252.0000       246.0000      1.000000    250.0000   0.000000      0
#> 25%-tile:   252.0000       249.0000      2.000000    251.0000   0.000000      0
#> Median:     253.0000       249.0000      2.000000    251.0000   1.000000      0
#> 75%-tile:   253.0000       250.0000      2.000000    251.0000   5.000000      0
#> 97.5%-tile: 254.0000       251.0000      4.000000    253.0000  26.000000      0
#> Maximum:    270.0000       255.0000     22.000000    256.0000 120.000000      0
#> Mean:       252.7575       249.1501      2.005361    251.1555   5.162474      0
#>             Expected_Errors
#> Minimum:         1.00000000
#> 2.5%-tile:       1.00000000
#> 25%-tile:        1.00000000
#> Median:          1.00000000
#> 75%-tile:        1.00000000
#> 97.5%-tile:      1.00000000
#> Maximum:         4.00000000
#> Mean:            0.07385095
#>               Length Overlap_Length Overlap_Start Overlap_End MisMatches Num_Ns
#> Minimum:    250.0000       232.0000      0.000000    248.0000   0.000000      0
#> 2.5%-tile:  252.0000       246.0000      1.000000    250.0000   0.000000      0
#> 25%-tile:   252.0000       249.0000      2.000000    251.0000   0.000000      0
#> Median:     253.0000       249.0000      2.000000    251.0000   1.000000      0
#> 75%-tile:   253.0000       250.0000      2.000000    251.0000   5.000000      0
#> 97.5%-tile: 254.0000       251.0000      4.000000    253.0000  26.000000      0
#> Maximum:    270.0000       255.0000     22.000000    256.0000 120.000000      0
#> Mean:       252.7575       249.1501      2.005361    251.1555   5.162474      0
#>             Expected_Errors
#> Minimum:         1.00000000
#> 2.5%-tile:       1.00000000
#> 25%-tile:        1.00000000
#> Median:          1.00000000
#> 75%-tile:        1.00000000
#> 97.5%-tile:      1.00000000
#> Maximum:         4.00000000
#> Mean:            0.07385095
#>               Length Overlap_Length Overlap_Start Overlap_End MisMatches Num_Ns
#> Minimum:    250.0000       232.0000      0.000000    248.0000   0.000000      0
#> 2.5%-tile:  252.0000       246.0000      1.000000    250.0000   0.000000      0
#> 25%-tile:   252.0000       249.0000      2.000000    251.0000   0.000000      0
#> Median:     253.0000       249.0000      2.000000    251.0000   1.000000      0
#> 75%-tile:   253.0000       250.0000      2.000000    251.0000   5.000000      0
#> 97.5%-tile: 254.0000       251.0000      4.000000    253.0000  26.000000      0
#> Maximum:    270.0000       255.0000     22.000000    256.0000 120.000000      0
#> Mean:       252.7575       249.1501      2.005361    251.1555   5.162474      0
#>             Expected_Errors
#> Minimum:         1.00000000
#> 2.5%-tile:       1.00000000
#> 25%-tile:        1.00000000
#> Median:          1.00000000
#> 75%-tile:        1.00000000
#> 97.5%-tile:      1.00000000
#> Maximum:         4.00000000
#> Mean:            0.07385095

# remove sample 'F3D0' to produce a scrap report
xdev_remove_samples(data = miseq, samples = c("F3D0"))
#> miseq_sop:
#> 
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  375    249      0        3     0      1.00
#> 2.5%-tile:       1  375    252      0        4     0   2694.30
#> 25%-tile:        1  375    252      0        4     0  26943.00
#> Median:          1  375    253      0        4     0  53886.00
#> 75%-tile:        1  375    253      0        5     0  80829.00
#> 97.5%-tile:      1  375    254      0        6     0 105077.70
#> Maximum:         1  375    256      0        6     0 107772.00
#> Mean:            1  375    252      0        4     0  53886.14
#> 
#> scrap_summary:
#>        type     trash_code unique total
#> 1  sequence remove_samples    101   109
#> 2       otu remove_samples     14    14
#> 3       asv remove_samples    101   109
#> 4 phylotype remove_samples      2     2
#> 
#> Number of unique seqs: 2324 
#> Total number of seqs: 107772 
#> 
#> Total number of samples: 18 
#> Total number of treatments: 2 
#> Total number of otus: 517 
#> Total number of otu bin classifications: 517 
#> Total number of asvs: 2324 
#> Total number of asv bin classifications: 2324 
#> Total number of phylotypes: 61 
#> Total number of phylotype bin classifications: 61 
#> Total number of sequence classifications: 2324 
#> Total number of resource references: 2 
#> Total number of custom reports: 2 
#> 

# summarize scrapped data -
# sequences and bins scrapped by removing the sample "F3D0"
miseq$summary(type = "scrap")
#>        type     trash_code unique total
#> 1  sequence remove_samples    101   109
#> 2       otu remove_samples     14    14
#> 3       asv remove_samples    101   109
#> 4 phylotype remove_samples      2     2
#>        type     trash_code unique total
#> 1  sequence remove_samples    101   109
#> 2       otu remove_samples     14    14
#> 3       asv remove_samples    101   109
#> 4 phylotype remove_samples      2     2
#>        type     trash_code unique total
#> 1  sequence remove_samples    101   109
#> 2       otu remove_samples     14    14
#> 3       asv remove_samples    101   109
#> 4 phylotype remove_samples      2     2
```
