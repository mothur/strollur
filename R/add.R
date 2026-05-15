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
#' 'sequence', 'resource_reference' 'metadata' and 'report'.
#'
#' @param report_type, a string containing the type of report you are adding.
#' Options include: 'metadata' and custom reports.
#'
#' @param table_names, named list used to indicate the names of the columns in
#' the table. By default:
#'
#' table_names <- list(sequence_name = "sequence_name",
#'                     comment = "comment",
#'                     sequence = "sequence",
#'                     reference_name = "name",
#'                     reference_vendor = "vendor",
#'                     reference_version = "version",
#'                     reference_usage = "usage",
#'                     reference_note = "note",
#'                     reference_documentation_url = "documentation_url",
#'                     reference_method_url = "method_url",
#'                     reference_parameter = "parameter",
#'                     reference_citation = "citation")
#'
#' In table_names, 'sequence_name' is a string containing the name of the column
#' in 'table' that contains the sequence names. It is used when you are adding
#' FASTA data. Default column name is 'sequence_name'.
#'
#' In table_names, 'sequence' is a string containing the name of the column in
#' 'table' that contains the sequence nucleotide strings. It is used when you
#' are adding FASTA data. Default column name is 'sequence'.
#'
#' In table_names, 'comment' is a string containing the name of the column in
#' 'table' that contains the sequence comments. It is used when you are adding
#' FASTA data. Default column name is 'comment'.
#'
#' In table_names, 'reference_vendor' is a string containing the name of the
#' column in 'table' that contains the reference vendor names. It is used when
#' you are adding reference data. Default column name is 'vendor'.
#
#' In table_names, 'reference_name' is a string containing the name of the
#' column in 'table' that contains the reference names. It is used when you are
#' adding reference data. Default column name is 'name'.
#'
#' In table_names, 'reference_version' is a string containing the name of the
#' column in 'table' that contains the reference versions. Default column name
#' is 'version'.
#'
#' In table_names, 'reference_usage' is a string containing the name of the
#' column in 'table' that contains the reference usages. Default column name is
#' 'usage'.
#'
#' In table_names, 'reference_note' is a string containing the name of the
#' column in 'table' that contains the reference notes. Default column name is
#'  'note'.
#'
#' In table_names, 'reference_method_url' is a string containing the name of the
#' column in 'table' that contains the reference method urls. Default column
#' name is 'method_url'.
#'
#' In table_names, 'reference_documentation_url' is a string containing the name
#' of the column in 'table' that contains the reference urls. Default column
#' name is 'documentation_url'.
#'
#' In table_names, 'reference_parameter' is a string containing the name of the
#' column in 'table' that contains the reference parameters. Default column name
#' is 'parameter'.
#'
#' In table_names, 'reference_citation' is a string containing the name of the
#' column in 'table' that contains the reference citations. Default column name
#' is 'citation'.
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
#' add(data = data, table = fasta_data, type = "sequence")
#'
#' # To add FASTA data with a resource reference
#'
#' # Create a new empty strollur object named 'example_dataset'
#' data <- new_dataset(dataset_name = "example_dataset")
#'
#' # Create a resource reference for the FASTA data silva_resource <-
#' silva_resource <- new_reference(
#'   vendor = "SILVA", name =
#'     "silva.bacteria.fasta", version = "1.38.1",
#'   usage = "alignment of sequences",
#'   note = "reference trimmed to V4 region", method_url =
#'     "https://mothur.org/blog/2024/SILVA-v138_2-reference-files/",
#'   documentation_url = "https://mothur.org/wiki/silva_reference_files/"
#' )
#'
#' # Add FASTA data with a resource reference
#'
#' add(
#'   data,
#'   table = fasta_data,
#'   type = "sequence",
#'   reference = silva_resource
#' )
#'
#' # Add contigs assembly report with a 'sequence_name' column named 'Name'
#'
#' contigs_report <- readRDS(strollur_example("miseq_contigs_report.rds"))
#'
#' add(
#'   data,
#'   table = contigs_report, type = "report",
#'   report_type = "contigs_report", list(sequence_name = "Name")
#' )
#'
#' # To add metadata related to your study
#'
#' metadata <- readRDS(strollur_example("miseq_metadata.rds"))
#'
#' add(data, table = metadata, type = "metadata")
#'
#' @return an updated
#'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#' @export
add <- function(data, table,
                type = "sequence",
                report_type = NULL,
                table_names = list(
                  sequence_name = "sequence_name",
                  sequence = "sequence",
                  comment = "comment",
                  reference_vendor = "vendor",
                  reference_name = "name",
                  reference_version = "version",
                  reference_usage = "usage",
                  reference_note = "note",
                  reference_method_url = "method_url",
                  reference_documentation_url = "documentation_url",
                  reference_parameter = "parameter",
                  reference_citation = "citation"
                ),
                reference = NULL,
                verbose = TRUE) {
  if (!inherits(data, "strollur")) {
    stop("data must be a strollur object.")
  }

  type <- as.character(substitute(type))

  default_tn <- list(
    sequence_name = "sequence_name",
    sequence = "sequence",
    comment = "comment",
    reference_vendor = "vendor",
    reference_name = "name",
    reference_version = "version",
    reference_usage = "usage",
    reference_note = "note",
    reference_method_url = "method_url",
    reference_documentation_url = "documentation_url",
    reference_parameter = "parameter",
    reference_citation = "citation"
  )

  table_names <- modifyList(default_tn, table_names)

  if (type == "sequence") {
    xdev_add_sequences(
      data = data, table = table,
      sequence_name = table_names[["sequence_name"]],
      sequence = table_names[["sequence"]],
      comment = table_names[["comment"]],
      reference = reference,
      verbose = verbose
    )
  } else if (type == "report") {
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
    xdev_add_report(
      data = data, table = table,
      type = type,
      verbose = verbose
    )
  } else if (type == "resource_reference") {
    xdev_add_references(
      data = data, table = table,
      name = table_names[["reference_name"]],
      vendor = table_names[["reference_vendor"]],
      version = table_names[["reference_version"]],
      usage = table_names[["reference_usage"]],
      note = table_names[["reference_note"]],
      method_url = table_names[["reference_method_url"]],
      documentation_url = table_names[["reference_documentation_url"]],
      parameter = table_names[["reference_parameter"]],
      citation = table_names[["reference_citation"]],
      verbose = verbose
    )
  } else {
    message <- paste0(
      type, " is not a valid 'type' for the add()",
      " function."
    )
    cli::cli_abort(message)
  }
  invisible(data)
}
