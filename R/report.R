#' @title Get a data.frame containing the given report in a
#' \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#' @name report
#' @rdname report
#' @description Get a data.frame containing the report or the phylo tree
#' requested. Reports include FASTA format, FASTQ format, sequences reports,
#' sequence_bin_assignments, sequence_taxonomy, sequence_tree, bin_taxonomy,
#' bin_representatives, sample_assignments, sample_tree, references, custom
#' reports, sequence_scrap, and bin_scrap in a
#' \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object.
#'
#' @param data a
#'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#'
#' @param type string containing the type of report you would like. Options
#' include: `fasta`, `fastq`, `sequence`, `sequence_bin_assignment`,
#' `sequence_taxonomy`, `sequence_tree`, `bin_taxonomy`, `bin_representative`,
#'  `sample_assignment`, `sample_tree`, `resource_reference`, `sequence_scrap`,
#' `bin_scrap`. If you have added custom reports for alignment,
#' contigs_assembly or chimeras, you can get those as well.
#'  Default = `sequence`.
#'
#' @param bin_type string containing the bin type you would like a bin_taxonomy
#' report for. Default = `otu`.
#'
#' @examples
#'
#' miseq <- strollur::miseq_sop_example()
#'
#' # To get the FASTA data
#'
#' strollur::report(data = miseq, type = "fasta") |> head(n = 5)
#'
#' # To get a report about the FASTA data
#'
#' strollur::report(data = miseq, type = "sequence") |> head(n = 5)
#'
#' # To get the sequence bin assignments
#'
#' strollur::report(data = miseq,
#'                  type = "sequence_bin_assignment", bin_type = "otu") |>
#'   head(n = 5)
#'
#' # To get the sample treatment assignments
#'
#' strollur::report(data = miseq, type = "sample_assignment")
#'
#' # To get a report about sequence classifications
#'
#' strollur::report(data = miseq, type = "sequence_taxonomy") |> head(n = 10)
#'
#' # To get a report about bin classifications for 'otu' data
#'
#' strollur::report(data = miseq,
#'                  type = "bin_taxonomy", bin_type = "otu") |> head(n = 10)
#'
#' # To get the 'otu' bin representative sequences
#'
#' strollur::report(
#'   data = miseq, type = "bin_representative",
#'   bin_type = "otu"
#' ) |> head(n = 5)
#'
#' # To get a report about the sequences removed during your analysis:
#'
#' strollur::report(data = miseq, type = "sequence_scrap")
#'
#' # To get a report about the "otu" bins removed during your analysis:
#'
#' strollur::report(data = miseq, type = "bin_scrap", bin_type = "otu")
#'
#' # To get the metadata associated with your data:
#'
#' metadata <- strollur::report(data = miseq, type = "metadata")
#'
#' # To get the resource references associated with your data:
#'
#' references <- strollur::report(data = miseq, type = "resource_reference")
#'
#' # To get our custom report containing the contigs assembly data:
#'
#' strollur::report(data = miseq, type = "contigs_report") |> head(n = 10)
#'
#' # To get the tree that relates your sequences:
#'
#' sequence_tree <- strollur::report(data = miseq, type = "sequence_tree")
#'
#' # To get the tree that relates your samples:
#'
#' sample_tree <- strollur::report(data = miseq, type = "sample_tree")
#'
#' @return data.frame or tree
#' @export
report <- function(data, type = "sequence", bin_type = "otu") {
  if (!inherits(data, "strollur")) {
    stop("data must be a strollur object.")
  }

  # allow for type and bin_type to be entered without ""
  type <- as.character(substitute(type))
  bin_type <- as.character(substitute(bin_type))

  if (type == "sequence_tree") {
    if (!is.null(data$sequence_tree)) {
      # prune tree if needed
      # seqs in tree and not in dataset
      extra_seqs <- setdiff(
        data$sequence_tree$tip.label,
        names(data, type = "sequence")
      )

      if (length(extra_seqs) != 0) {
        # if tree contains "extra" names, prune the tree
        data$sequence_tree <- drop.tip(data$sequence_tree,
          tip = extra_seqs
        )
      }
    }
    return(data$sequence_tree)
  } else if (type == "sample_tree") {
    if (!is.null(data$sample_tree)) {
      # prune tree if needed
      # samples in tree and not in dataset
      extra_samples <- setdiff(
        data$sample_tree$tip.label,
        names(data, type = "sample")
      )

      if (length(extra_samples) != 0) {
        # if tree contains "extra" samples, prune the tree
        data$sample_tree <- drop.tip(data$sample_tree,
          tip = extra_samples
        )
      }
    }
    return(data$sample_tree)
  }
  xdev_report(data, type, bin_type)
}
