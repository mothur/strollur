#' @title Add sequences, reports, metadata or resource references to a
#'   \link{strollur} object
#' @name add
#' @rdname add
#' @description
#' Add sequences, reports, metadata or resource references to a
#' \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#'
#' @param data, a
#'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#'
#' @param table, a data.frame containing the data you wish to add.
#'
#' @param type, a string containing the type of data. Options include:
#' 'sequences', 'references' 'metadata' and 'reports'.
#'
#' @param report_type, a string containing the type of report you are adding.
#' Options include: 'metadata' and custom reports.
#'
#' @param table_names, named list used to indicate the names of the columns in
#' the table. By default:
#'
#' table_names <- list(sequence_name = "sequence_names",
#'                     comment = "comments",
#'                     sequence = "sequences",
#'                     reference_name = "reference_names",
#'                     reference_version = "reference_versions",
#'                     reference_usage = "reference_usages",
#'                     reference_note = "reference_notes",
#'                     reference_url = "reference_urls")
#'
#' In table_names, 'sequence_name' is a string containing the name of the column
#' in 'table' that contains the sequence names. It is used when you are adding
#' FASTA data. Default column name is 'sequence_names'.
#'
#' In table_names, 'sequence' is a string containing the name of the column in
#' 'table' that contains the sequence nucleotide strings. It is used when you
#' are adding FASTA data. Default column name is 'sequences'.
#'
#' In table_names, 'comment' is a string containing the name of the column in
#' 'table' that contains the sequence comments. It is used when you are adding
#' FASTA data. Default column name is 'comments'.
#'
#' In table_names, 'reference_name' is a string containing the name of the
#' column in 'table' that contains the reference names. It is used when you are
#' adding reference data. Default column name is 'reference_names'.
#'
#' In table_names, 'reference_version' is a string containing the name of the
#' column in 'table' that contains the reference versions. Default column name
#' is 'reference_versions'.
#'
#' In table_names, 'reference_usage' is a string containing the name of the
#' column in 'table' that contains the reference usages. Default column name is
#' 'reference_usages'.
#'
#' In table_names, 'reference_note' is a string containing the name of the
#' column in 'table' that contains the reference notes. Default column name is
#'  'reference_notes'.
#'
#' In table_names, 'reference_url' is a string containing the name of the column
#' in 'table' that contains the reference urls. Default column name is
#'  'reference_urls'.
#'
#' @param reference, a list created by the function [new_reference]. Optional.
#'
#' @param verbose, boolean indicating whether or not you want progress messages.
#' Default = TRUE.
#'
#' @examples
#'
#' # Create a new empty strollur object named 'example_dataset'
#' data <- new_dataset(dataset_name = "example_dataset")
#'
#' # Read FASTA data into data.frame
#' fasta_data <- read_fasta(fasta = strollur_example("final.fasta.gz"))
#'
#' # Add FASTA sequence data
#' add(data = data, table = fasta_data, type = "sequences")
#'
#' # To add FASTA data with a resource reference
#'
#' # Create a new empty strollur object named 'example_dataset'
#' data <- new_dataset(dataset_name = "example_dataset")
#'
#' # Create a resource reference for the FASTA data
#' resource_url <- "https://mothur.org/wiki/silva_reference_files/"
#' resource_reference <-
#'   new_reference(
#'     reference_name = "silva.bacteria.fasta",
#'     reference_version = "1.38.1",
#'     reference_usage = "alignment by mothur2 v1.0",
#'     reference_note = "default options",
#'     reference_url = resource_url
#'   )
#'
#' # Add FASTA data with a resource reference
#'
#' add(
#'   data,
#'   table = fasta_data,
#'   type = "sequences",
#'   reference = resource_reference
#' )
#'
#' # Add contigs assembly report with a 'sequence_name' column named 'Name'
#'
#' contigs_report <- readRDS(strollur_example("miseq_contigs_report.rds"))
#'
#' add(
#'   data,
#'   table = contigs_report, type = "reports",
#'   report_type = "contigs_report", list(sequence_name = "Name")
#' )
#'
#' # To add metadata related to your study
#'
#' metadata <- readRDS(strollur_example("miseq_metadata.rds"))
#'
#' add(data, table = metadata, type = "metadata")
#'
#' @return double - The number of items added
#' @export
add <- function(data, table,
                type = "sequences",
                report_type = NULL,
                table_names = list(
                  sequence_name = "sequence_names",
                  sequence = "sequences",
                  comment = "comments",
                  reference_name = "reference_names",
                  reference_version = "reference_versions",
                  reference_usage = "reference_usages",
                  reference_note = "reference_notes",
                  reference_url = "reference_urls"
                ),
                reference = NULL,
                verbose = TRUE) {
  if (!inherits(data, "strollur")) {
    stop("data must be a strollur object.")
  }

  type <- as.character(substitute(type))

  default_tn <- list(
    sequence_name = "sequence_names",
    sequence = "sequences",
    comment = "comments",
    reference_name = "reference_names",
    reference_version = "reference_versions",
    reference_usage = "reference_usages",
    reference_note = "reference_notes",
    reference_url = "reference_urls"
  )

  table_names <- modifyList(default_tn, table_names)

  num_added <- 0
  if (type == "sequences") {
    num_added <- xdev_add_sequences(
      data = data, table = table,
      sequence_name = table_names[["sequence_name"]],
      sequence = table_names[["sequence"]],
      comment = table_names[["comment"]],
      reference = reference,
      verbose = verbose
    )
  } else if (type == "reports") {
    num_added <- 1

    if (!is.null(report_type)) {
      xdev_add_report(
        data = data, table = table,
        type = report_type,
        sequence_name = table_names[["sequence_name"]],
        verbose
      )
    } else {
      cli::cli_abort("'report_type' is required when adding a report.")
    }
  } else if (type == "metadata") {
    num_added <- 1
    xdev_add_report(
      data = data, table = table,
      type = type,
      verbose = verbose
    )
  } else if (type == "references") {
    num_added <- xdev_add_references(
      data = data, table = table,
      reference_name = table_names[["reference_name"]],
      reference_version = table_names[["reference_version"]],
      reference_usage = table_names[["reference_usage"]],
      reference_note = table_names[["reference_note"]],
      reference_url = table_names[["reference_url"]],
      verbose = verbose
    )
  } else {
    message <- paste0(
      type, " is not a valid 'type' for the add()",
      " function."
    )
    cli::cli_abort(message)
  }
  num_added
}
