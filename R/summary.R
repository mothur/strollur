# Define an S3 generic - this allows for the additional parameters to summary
summary <- function(x, type = "sequences",
                    bin_type = "otu",
                    verbose = TRUE, ...) {
  UseMethod("summary", x)
}

#' @title Summarize the sequences data, custom reports, and scrapped data in a
#' \link{strollur::strollur} object
#' @name summary
#' @rdname summary
#' @description
#' Summarize the sequences data, custom reports, and scrapped data in a
#' \link{strollur::strollur} object
#'
#' @param data, a \link{strollur::strollur} object
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
  if ("strollur" %in% class(data)) {
    dataset_summary <- xdev_summarize(data, type, report_type)
    if (verbose) {
      print(dataset_summary)
    }
    dataset_summary
  } else {
    base::summary(data)
  }
}
