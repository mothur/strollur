# Define an S3 generic - this allows for the additional parameters to summary
summary <- function(x, type = "sequence",
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
#' Options include: "sequence", "report" and "scrap". Default = "sequence".
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
#' summary(data = miseq, type = "sequence")
#'
#' # summarize contigs_report
#' summary(data = miseq, type = "report", report_type = "contigs_report")
#'
#' # remove sample 'F3D0' to produce a scrap report
#' xdev_remove_samples(data = miseq, samples = c("F3D0"))
#'
#' # summarize FASTA data after removal of sample F3D0
#' summary(data = miseq, type = "sequence")
#'
#' # summarize scrapped data -
#' # sequences and bins scrapped by removing the sample "F3D0"
#' summary(data = miseq, type = "scrap")
#'
#' @return data.frame
#' @export
summary <- function(data, type = "sequence",
                    report_type = NULL, verbose = TRUE) {
  if (!inherits(data, "strollur")) {
    stop("data must be a strollur object.")
  }

  type <- as.character(substitute(type))
  if (!is.null(report_type)) {
    report_type <- as.character(substitute(report_type))
  }

  if (all(type != c("sequence", "report", "scrap"))) {
    stop(paste0(
      type, " is not a valid type option. Options include: ",
      "'sequence', 'report' and 'scrap'."
    ))
  }

  if (!is.null(report_type)) {
    if (all(report_type != xdev_names(data = data, type = "report"))) {
      stop(paste0(report_type), " is not a valid report_type option.")
    }
  }

  dataset_summary <- ""
  if (type == "sequence") {
    dataset_summary <- generate_sequence_report(data)
  } else if (type == "report") {
    dataset_summary <- generate_report(data, report_type)
  } else { # scrap reports still use cpp
    dataset_summary <- xint_get_scrap_summary(data)
  }

  if (verbose) {
    print(dataset_summary)
  }
  dataset_summary
}

generate_sequence_report <- function(dataset) {
  report <- report(dataset)
  abunds <- abundance(dataset)
  report <- cbind(report, abundance = abunds)

  desired_quantiles <- c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1)
  desired_tags <- c(
    "Minimum:", "2.5%-tile:", "25%-tile:", "Median:",
    "75%-tile:", "97.5%-tile:", "Maximum:", "Mean:"
  )

  report_summary <- report |>
    reframe(
      stat = desired_tags,
      across(
        where(is.numeric),
        ~ c(
          quantile(.x, probs = desired_quantiles, na.rm = TRUE),
          mean(.x, na.rm = TRUE)
        ),
        .names = "{.col}"
      )
    )

  # remove weighted counts
  report_summary$abundance.abundance <- NULL

  total_seqs <- count(dataset, "sequence")
  num_seqs <- total_seqs * desired_quantiles
  num_seqs[[1]] <- 1 # minimum seqs should be 1
  num_seqs <- c(num_seqs, mean(num_seqs))
  report_summary <- cbind(report_summary, num_seqs)
  rownames(report_summary) <- report_summary$stat
  colnames(report_summary) <- c(
    "stat", "starts", "ends", "nbases", "ambigs",
    "polymers", "numns", "numseqs"
  )
  report_summary[, -1]
}

generate_report <- function(dataset, report_type) {
  desired_quantiles <- c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1)
  desired_tags <- c(
    "Minimum:", "2.5%-tile:", "25%-tile:", "Median:",
    "75%-tile:", "97.5%-tile:", "Maximum:", "Mean:"
  )

  result <- xdev_report(dataset, report_type) |>
    reframe(
      stat = desired_tags,
      across(
        where(is.numeric),
        ~ c(
          ceiling(quantile(.x, probs = desired_quantiles, na.rm = TRUE)),
          mean(.x, na.rm = TRUE)
        ),
        .names = "{.col}"
      )
    )
  rownames(result) <- result$stat
  result[, -1]
}
