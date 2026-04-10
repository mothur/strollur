#' @title Get a data.frame containing the given report in a
#' \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#' @name report
#' @rdname report
#' @description
#' Get a data.frame containing the report. Reports include FASTA format,
#' sequences reports, sequence_bin_assignments, sequence_taxonomy, bin_taxonomy,
#' bin_representatives, sample_assignments, metadata, references,
#' sequence_scrap, and bin_scrap in a
#' \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object.
#'
#' @param data, a
#'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#'
#' @param type, string containing the type of report you would like. Options
#' include: "fasta", "sequences", "sequence_bin_assignments",
#' "sequence_taxonomy", "bin_taxonomy", "bin_representatives",
#'  "sample_assignments", "metadata", "references", "sequence_scrap",
#' "bin_scrap". If you have added custom reports for alignment,
#' contigs_assembly or chimeras, you can get those as well.
#'  Default = "sequences".
#'
#' @param bin_type, string containing the bin type you would like a bin_taxonomy
#' report for. Default = "otu".
#'
#' @examples
#'
#' miseq <- miseq_sop_example()
#'
#' # To get the FASTA data
#'
#' report(data = miseq, type = "fasta") |> head(n = 5)
#'
#' # To get a report about the FASTA data
#'
#' report(data = miseq, type = "sequences") |> head(n = 5)
#'
#' # To get the sequence bin assignments
#'
#' report(data = miseq, type = "sequence_bin_assignments", bin_type = "otu") |>
#'   head(n = 5)
#'
#' # To get the sample treatment assignments
#'
#' report(data = miseq, type = "sample_assignments")
#'
#' # To get a report about sequence classifications
#'
#' report(data = miseq, type = "sequence_taxonomy") |> head(n = 10)
#'
#' # To get a report about bin classifications for 'otu' data
#'
#' report(data = miseq, type = "bin_taxonomy", bin_type = "otu") |> head(n = 10)
#'
#' # To get the 'otu' bin representative sequences
#'
#' report(
#'   data = miseq, type = "bin_representatives",
#'   bin_type = "otu"
#' ) |> head(n = 5)
#'
#' # To get a report about the sequences removed during your analysis:
#'
#' report(data = miseq, type = "sequence_scrap")
#'
#' # To get a report about the "otu" bins removed during your analysis:
#'
#' report(data = miseq, type = "bin_scrap", bin_type = "otu")
#'
#' # To get the metadata associated with your data:
#'
#' metadata <- report(data = miseq, type = "metadata")
#'
#' # To get the resource references associated with your data:
#'
#' references <- report(data = miseq, type = "references")
#'
#' # To get our custom report containing the contigs assembly data:
#'
#' report(data = miseq, type = "contigs_report") |> head(n = 10)
#'
#' @return data.frame
#' @export
report <- function(data, type = "sequences", bin_type = "otu") {
  if (!inherits(data, "strollur")) {
    stop("data must be a strollur object.")
  }

  # allow for type and bin_type to be entered without ""
  type <- as.character(substitute(type))
  bin_type <- as.character(substitute(bin_type))

  xdev_report(data, type, bin_type)
}
