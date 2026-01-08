#' @title read_qiime2_metadata
#' @description
#' Read a \href{https://qiime2.org}{qiime2} .tsv table containing metadata.
#'
#' @param metadata file name, a qiime2 .tsv file containing metadata about your
#' analysis.
#'
#' @examples
#'
#' metadata <- read_qiime2_metadata(rdataset_example(
#'   "sample_metadata.tsv"
#' ))
#'
#' @return A data.frame containing metadata
#' @export
read_qiime2_metadata <- function(metadata) {
  if (!file.exists(metadata)) {
    .abort_nonexistant_file(metadata)
  }

  # read entire file
  file_conn <- file(metadata)
  file_data <- suppressWarnings(readLines(file_conn))
  num_lines <- length(file_data)

  metadata <- data.frame()

  if (num_lines > 2) {
    # parse header line
    headers <- .split_at_char(file_data[[1]], "\t")
    headers[1] <- "SampleID"

    # look at second line and check for "#q2:types"
    ignore_line <- grepl("#q2:types", file_data[[2]])

    # create data.frame from rest of table
    if (ignore_line) {
      metadata <- read.table(file_conn,
        header = FALSE, col.names = headers,
        skip = 2, sep = "\t", check.names = FALSE
      )
    } else {
      metadata <- read.table(file_conn,
        header = FALSE, col.names = headers,
        skip = 1, sep = "\t", check.names = FALSE
      )
    }
  }

  metadata
}
