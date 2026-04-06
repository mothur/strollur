# Define an S3 generic - this allows for the additional parameters to summary
summary <- function(x, type = "sequences",
                    bin_type = "otu",
                    verbose = TRUE, ...) {
  UseMethod("summary", x)
}

#' @title Summarize the sequences data, custom reports, and scrapped data in a
#' \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#' @name summary
#' @rdname summary
#' @description
#' Summarize the sequences data, custom reports, and scrapped data in a
#' \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#'
#' @param data, a
#'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#'
#' @param type, string containing the type of data you want the number of.
#' Options include: "sequences", "reports" and "scrap". Default = "sequences".
#'
#' @param report_type, string containing the report type you would summarized.
#' For example, the miseq_sop_example includes contigs assembly data and can be
#' accessed with report_type = "contigs_report". Default = NULL.
#'
#' @param verbose, boolean indicating whether or not you want progress messages.
#' Default = TRUE.
#'
#' @examples
#'
#' miseq <- miseq_sop_example()
#'
#' # To get the summary of your FASTA data
#' summary(data = miseq, type = "sequences")
#'
#' # summarize contigs_report
#' summary(data = miseq, type = "reports", report_type = "contigs_report")
#'
#' # remove sample 'F3D0' to produce a scrap report
#' xdev_remove_samples(data = miseq, samples = c("F3D0"))
#'
#' # summarize FASTA data after removal of sample F3D0
#' summary(data = miseq, type = "sequences")
#'
#' # summarize scrapped data -
#' # sequences and bins scrapped by removing the sample "F3D0"
#' summary(data = miseq, type = "scrap")
#'
#' @return data.frame
#' @export
summary <- function(data, type = "sequences",
                    report_type = NULL, verbose = TRUE) {
  if(type == "sequences") {
    return(generate_sequence_report(data))
  }
  else if(type == "report" && report_type == "contigs_report") {
    return(generate_contig_report(data))
  }
  # otherwise, we can summarize it differently
  if ("dataset" %in% class(data)) {
    dataset_summary <- xdev_summarize(data, type, report_type)
    if (verbose) {
      print(dataset_summary)
    }
    dataset_summary
  } else {
    base::summary(data)
  }
}

generate_sequence_report <- function(dataset) {
  dataset <- dat
  report <- report(dataset) 
  abunds <- abundance(dataset)
  report <- cbind(report, abundance = abunds)

  desired_quantiles <- c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1)
  desired_tags <- c("Minimum:", "2.5%-tile:", "25%-tile:", "Median:",
  "75%-tile:", "97.5%-tile:", "Maximum:", "Mean:")

  report_summary <- report |>
    reframe(stat = desired_tags,
      across(
      where(is.numeric),
      ~ c(quantile(.x, probs = desired_quantiles, na.rm = TRUE),
      mean(.x, na.rm = TRUE)),
      .names = "{.col}"
    )
  )

  # remove weighted counts
  report_summary$abundance.abundances <- NULL
  report_summary$abundance.sequence_names <- NULL


  total_seqs <- count(dataset, type)
  num_seqs <- total_seqs * desired_quantiles + 1 # minimum seqs should be 1
  num_seqs <- c(num_seqs, mean(num_seqs))
  report_summary <- cbind(report_summary, num_seqs)
  report_summary
}

generate_contig_report <- function(dataset) {
  desired_quantiles <- c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1)
  desired_tags <- c("Minimum:", "2.5%-tile:", "25%-tile:", "Median:",
  "75%-tile:", "97.5%-tile:", "Maximum:", "Mean:")

  report(dat, "contigs_report") |>
      reframe(stat = desired_tags,
        across(
        where(is.numeric),
        ~ c(quantile(.x, probs = desired_quantiles, na.rm = TRUE),
        mean(.x, na.rm = TRUE)),
        .names = "{.col}"
      )
    )
}
