#' @title Assign sequence abundances, sequence classifications, bins, bin
#' representative sequences, bin classifications or treatments to a
#' \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#' @name assign
#' @rdname assign
#' @description
#' Assign sequence abundances, sequence classifications, bins, bin
#' representative sequences, bin classifications or treatments to a
#' \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#'
#' @param data, a
#'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#'
#' @param table, a data.frame containing the data you wish to assign
#'
#' @param type, a string containing the type of data. Options include:
#' 'sequence_abundance', 'sequence_taxonomy', 'bins',
#'  'bin_representatives', 'bin_taxonomy' and 'treatments'.
#'  Default = "bins".
#'
#' @param bin_type, string containing the bin type you would like the number of
#' bins for. Default = "otu".
#'
#' @param table_names, named list used to indicate the names of the columns in
#' the table. By default:
#'
#' table_names <- list(sequence_name = "sequence_names",
#'                     abundance = "abundances",
#'                     sample = "samples",
#'                     treatment = "treatments",
#'                     taxonomy = "taxonomies",
#'                     bin_name = "bin_names")
#'
#' In table_names, 'sequence_name' is a string containing the name of the column
#' in 'table' that contains the sequence names. Default column name is
#' 'sequence_names'.
#'
#' In table_names, 'abundance' is a string containing the name of the column in
#' 'table' that contains the abundances. Default column name is 'abundances'.
#'
#' In table_names, 'sample' is a string containing the name of the column in
#' 'table' that contains the samples. Default column name is 'samples'.
#'
#' In table_names, 'treatment' is a string containing the name of the
#' column in 'table' that contains the treatment names. Default column name is
#'  'treatments'.
#'
#' In table_names, 'taxonomy' is a string containing the name of the
#' column in 'table' that contains the classifications. Default column name
#' is 'taxonomies'.
#'
#' In table_names, 'bin_name' is a string containing the name of the
#' column in 'table' that contains the bin names. Default column name is
#' 'bin_names'.
#'
#' @param reference, a list created by the function [new_reference]. Optional.
#' @param verbose, boolean indicating whether or not you want progress messages.
#' Default = TRUE.
#'
#' @examples
#'
#' # Assign sequence classifications
#'
#' # create a new empty strollur object named 'example_dataset'
#' data <- new_dataset(dataset_name = "example_dataset")
#'
#' sequence_classifications <- read_mothur_taxonomy(strollur_example(
#'   "final.taxonomy.gz"
#' ))
#'
#' assign(
#'   data,
#'   table = sequence_classifications, type = "sequence_taxonomy"
#' )
#'
#' # Assigning bins
#'
#' # read mothur's otu list file into data.frame
#' otu_data <- read_mothur_list(list = strollur_example(
#'   "final.opti_mcc.list.gz"
#' ))
#'
#' # read mothur's asv list file into data.frame
#' asv_data <- read_mothur_list(list = strollur_example(
#'   "final.asv.list.gz"
#' ))
#'
#' # read mothur's phylotype list file into data.frame
#' phylo_data <- read_mothur_list(list = strollur_example(
#'   "final.tx.list.gz"
#' ))
#'
#' # read otu bin representative sequences into a data.frame
#' bin_reps <- readRDS(strollur_example("miseq_representative_sequences.rds"))
#'
#' # assign 'otu' bins using sequence names
#' assign(data, table = otu_data, bin_type = "otu")
#'
#' # assign 'asv' bins using sequence names
#' assign(data, table = asv_data, bin_type = "asv")
#'
#' # assign 'phylotype' bins using sequence names
#' assign(data, table = phylo_data, bin_type = "phylotype")
#'
#' # assign 'otu' bin representative sequences
#' assign(data, table = bin_reps, type = "bin_representatives")
#'
#' # To assign abundance only bins
#'
#' # create a new empty strollur object named 'example_dataset'
#' data <- new_dataset(dataset_name = "example_dataset")
#'
#' # read mothur's shared file
#' otu_data <- read_mothur_shared(strollur_example("final.opti_mcc.shared"))
#'
#' # assign abundance only otus parsed by sample
#' assign(data, table = otu_data, bin_type = "otu")
#'
#' # Assigning bin classifications
#'
#' # read bin taxonomies
#' otu_data <- read_mothur_cons_taxonomy(strollur_example(
#'   "final.cons.taxonomy"
#' ))
#'
#' # assign otu consensus taxonomies
#' assign(
#'   data,
#'   table = otu_data,
#'   type = "bin_taxonomy", bin_type = "otu"
#' )
#'
#' # Assign treatments
#'
#' sample_assignments <- readRDS(strollur_example("miseq_sample_design.rds"))
#'
#' assign(data, table = sample_assignments, type = "treatments")
#'
#' @return double - The number of items assigned
#' @export
assign <- function(data, table,
                   type = "bins",
                   bin_type = "otu",
                   table_names = list(
                     sequence_name = "sequence_names",
                     abundance = "abundances",
                     sample = "samples",
                     treatment = "treatments",
                     taxonomy = "taxonomies",
                     bin_name = "bin_names"
                   ),
                   reference = NULL,
                   verbose = TRUE) {
  if (!inherits(data, "strollur")) {
    stop("data must be a strollur object.")
  }

  default_tn <- list(
    sequence_name = "sequence_names",
    abundance = "abundances",
    sample = "samples",
    treatment = "treatments",
    taxonomy = "taxonomies",
    bin_name = "bin_names"
  )

  table_names <- modifyList(default_tn, table_names)

  # allow for type and bin_type to be entered without ""
  type <- as.character(substitute(type))
  bin_type <- as.character(substitute(bin_type))

  num <- 0
  if (type == "bins") {
    num <- xdev_assign_bins(
      data = data, table = table,
      bin_type = bin_type,
      reference = reference,
      bin_name = table_names[["bin_name"]],
      abundance = table_names[["abundance"]],
      sample = table_names[["sample"]],
      sequence_name = table_names[["sequence_name"]],
      verbose = verbose
    )
  } else if (type == "bin_taxonomy") {
    num <- xdev_assign_bin_taxonomy(
      data = data, table = table,
      bin_type = bin_type,
      reference = reference,
      bin_name = table_names[["bin_name"]],
      taxonomy = table_names[["taxonomy"]],
      verbose = verbose
    )
  } else if (type == "bin_representatives") {
    num <- xdev_assign_bin_representative_sequences(
      data = data, table = table,
      bin_type = bin_type,
      reference = reference,
      bin_name = table_names[["bin_name"]],
      sequence_name = table_names[["sequence_name"]],
      verbose = verbose
    )
  } else if (type == "sequence_taxonomy") {
    num <- xdev_assign_sequence_taxonomy(
      data = data, table = table,
      reference = reference,
      sequence_name = table_names[["sequence_name"]],
      taxonomy = table_names[["taxonomy"]],
      verbose = verbose
    )
  } else if (type == "treatments") {
    num <- xdev_assign_treatments(
      data = data, table = table,
      sample = table_names[["sample"]],
      treatment = table_names[["treatment"]],
      verbose = verbose
    )
  } else if (type == "sequence_abundance") {
    num <- xdev_assign_sequence_abundance(
      data = data, table = table,
      sequence_name = table_names[["sequence_name"]],
      abundance = table_names[["abundance"]],
      sample = table_names[["sample"]],
      treatment = table_names[["treatment"]],
      verbose = verbose
    )
  } else {
    message <- paste0(
      type, " is not a valid 'type' for the assign()",
      " function."
    )
    cli::cli_abort(message)
  }
  num
}
