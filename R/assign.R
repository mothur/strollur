#' @title Assign sequence abundances, sequence classifications, bins, bin
#' representative sequences, bin classifications or treatments to a
#' \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#' @name assign
#' @rdname assign
#'
#' @description
#' Assign sequence abundances, sequence classifications, bins, bin
#' representative sequences, bin classifications or treatments to a
#' \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#'
#' @param data a
#'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#'
#' @param table a data.frame containing the data you wish to assign
#'
#' @param type a string containing the type of data. Options include:
#' `sequence_abundance`, `sequence_taxonomy`, `bin`, `quality`,
#'  `bin_representative`, `bin_taxonomy` and `treatment`.
#'  Default = `bin`.
#'
#' @param bin_type string containing the bin type you would like the number of
#' bins for. Default = `otu`.
#'
#' @param table_names named list used to indicate the names of the columns in
#' the table. By default:
#'
#' table_names <- list(sequence_name = "sequence_name",
#'                     quality_score = "quality_score",
#'                     abundance = "abundance",
#'                     sample = "sample",
#'                     treatment = "treatment",
#'                     taxonomy = "taxonomy",
#'                     level = "level",
#'                     confidence = "confidence",
#'                     bin_name = "bin_name")
#'
#' In table_names, `sequence_name` is a string containing the name of the column
#' in 'table' that contains the sequence names. Default column name is
#' 'sequence_name'.
#'
#' In table_names, `quality_score` is a string containing the name of the column
#' in 'table' that contains the sequence quality scores. It is used when you are
#' adding quality data. Default column name is 'quality_score'.
#'
#' In table_names, `abundance` is a string containing the name of the column in
#' 'table' that contains the abundances. Default column name is 'abundance'.
#'
#' In table_names, `sample` is a string containing the name of the column in
#' 'table' that contains the samples. Default column name is 'sample'.
#'
#' In table_names, `treatment` is a string containing the name of the
#' column in 'table' that contains the treatment names. Default column name is
#'  'treatment'.
#'
#' In table_names, `taxonomy` is a string containing the name of the
#' column in 'table' that contains the classifications. Default column name
#' is 'taxonomy'.
#'
#' In table_names, `level` is a string containing the name of the column
#' in the data.frame that contains the classifications. Default column name is
#' `level`.
#'
#' In table_names, `confidence` is a string containing the name of the column in
#' the data.frame that contains the classifications confidence scores. Default
#' column name is `confidence`.
#'
#' In table_names, `bin_name` is a string containing the name of the
#' column in 'table' that contains the bin names. Default column name is
#' 'bin_name'.
#'
#' @param reference a list created by the function [new_reference]. Optional.
#' @param verbose boolean indicating whether or not you want progress messages.
#' Default = TRUE.
#'
#' @seealso [taxonomy_table_format]
#'
#' @examples
#'
#' # Assign sequence classifications
#'
#' # create a new empty strollur object named 'example_dataset'
#' data <- strollur::new_dataset(dataset_name = "example_dataset")
#'
#' sequence_classifications <- strollur::read_mothur_taxonomy(strollur_example(
#'   "final.taxonomy.gz"
#' ))
#'
#' strollur::assign(
#'   data,
#'   table = sequence_classifications, type = "sequence_taxonomy"
#' )
#'
#' # Assigning bins
#'
#' # read mothur's otu list file into data.frame
#' otu_data <- strollur::read_mothur_list(list = strollur_example(
#'   "final.opti_mcc.list.gz"
#' ))
#'
#' # read mothur's asv list file into data.frame
#' asv_data <- strollur::read_mothur_list(list = strollur_example(
#'   "final.asv.list.gz"
#' ))
#'
#' # read mothur's phylotype list file into data.frame
#' phylo_data <- strollur::read_mothur_list(list = strollur_example(
#'   "final.tx.list.gz"
#' ))
#'
#' # read otu bin representative sequences into a data.frame
#' bin_reps <- readRDS(strollur_example("miseq_representative_sequences.rds"))
#'
#' # assign 'otu' bins using sequence names
#' strollur::assign(data, table = otu_data, bin_type = "otu")
#'
#' # assign 'asv' bins using sequence names
#' strollur::assign(data, table = asv_data, bin_type = "asv")
#'
#' # assign 'phylotype' bins using sequence names
#' strollur::assign(data, table = phylo_data, bin_type = "phylotype")
#'
#' # assign 'otu' bin representative sequences
#' strollur::assign(data, table = bin_reps, type = "bin_representative")
#'
#' # To assign abundance only bins
#'
#' # create a new empty strollur object named 'example_dataset'
#' data <- strollur::new_dataset(dataset_name = "example_dataset")
#'
#' # read mothur's shared file
#' otu_data <-
#'   strollur::read_mothur_shared(strollur_example("final.opti_mcc.shared"))
#'
#' # assign abundance only otus parsed by sample
#' strollur::assign(data, table = otu_data, bin_type = "otu")
#'
#' # Assigning bin classifications
#'
#' # read bin taxonomies
#' otu_data <- strollur::read_mothur_cons_taxonomy(strollur_example(
#'   "final.cons.taxonomy"
#' ))
#'
#' # assign otu consensus taxonomies
#' strollur::assign(
#'   data,
#'   table = otu_data,
#'   type = "bin_taxonomy", bin_type = "otu"
#' )
#'
#' # Assign treatments
#'
#' sample_assignments <- readRDS(strollur_example("miseq_sample_design.rds"))
#'
#' strollur::assign(data, table = sample_assignments, type = "treatment")
#'
#' @return an updated
#'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#' @export
assign <- function(data, table,
                   type = "bin",
                   bin_type = "otu",
                   table_names = list(
                     sequence_name = "sequence_name",
                     quality_score = "quality_score",
                     abundance = "abundance",
                     sample = "sample",
                     treatment = "treatment",
                     taxonomy = "taxonomy",
                     level = "level",
                     confidence = "confidence",
                     bin_name = "bin_name"
                   ),
                   reference = NULL,
                   verbose = TRUE) {
  if (!inherits(data, "strollur")) {
    stop("data must be a strollur object.")
  }

  default_tn <- list(
    sequence_name = "sequence_name",
    quality_score = "quality_score",
    abundance = "abundance",
    sample = "sample",
    treatment = "treatment",
    taxonomy = "taxonomy",
    level = "level",
    confidence = "confidence",
    bin_name = "bin_name"
  )

  table_names <- modifyList(default_tn, table_names)

  # allow for type and bin_type to be entered without ""
  type <- as.character(substitute(type))
  bin_type <- as.character(substitute(bin_type))

  num <- 0
  if (type == "bin") {
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
    # is this a tidy table?
    if (table_names[["level"]] %in% names(table)) {
      num <- xdev_assign_bin_taxonomy_tidy(
        data = data, table = table,
        bin_type = bin_type,
        reference = reference,
        bin_name = table_names[["bin_name"]],
        level = table_names[["level"]],
        taxonomy = table_names[["taxonomy"]],
        confidence = table_names[["confidence"]],
        verbose = verbose
      )
    } else {
      num <- xdev_assign_bin_taxonomy(
        data = data, table = table,
        bin_type = bin_type,
        reference = reference,
        bin_name = table_names[["bin_name"]],
        taxonomy = table_names[["taxonomy"]],
        verbose = verbose
      )
    }
  } else if (type == "bin_representative") {
    num <- xdev_assign_bin_representative_sequences(
      data = data, table = table,
      bin_type = bin_type,
      reference = reference,
      bin_name = table_names[["bin_name"]],
      sequence_name = table_names[["sequence_name"]],
      verbose = verbose
    )
  } else if (type == "sequence_taxonomy") {
    # is this a tidy table?
    if (table_names[["level"]] %in% names(table)) {
      num <- xdev_assign_sequence_taxonomy_tidy(
        data = data, table = table,
        reference = reference,
        sequence_name = table_names[["sequence_name"]],
        level = table_names[["level"]],
        taxonomy = table_names[["taxonomy"]],
        confidence = table_names[["confidence"]],
        verbose = verbose
      )
    } else {
      num <- xdev_assign_sequence_taxonomy(
        data = data, table = table,
        reference = reference,
        sequence_name = table_names[["sequence_name"]],
        taxonomy = table_names[["taxonomy"]],
        verbose = verbose
      )
    }
  } else if (type == "treatment") {
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
  } else if (type == "quality") {
    num <- xdev_assign_sequence_quality_scores(
      data = data, table = table,
      sequence_name = table_names[["sequence_name"]],
      quality_score = table_names[["quality_score"]],
      verbose = verbose
    )
  } else {
    message <- paste0(
      type, " is not a valid 'type' for the assign()",
      " function."
    )
    cli::cli_abort(message)
  }
  invisible(data)
}
