#' @title write_quality
#' @description
#'
#' Write a file containing sequence quality scores
#'
#' @param data a
#'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#' @param filename a string containing the name of the output file. Default =
#' 'dataset_name'.qual
#' @examples
#'
#' table <- strollur::read_fastq(strollur_example("tiny.fastq.gz"))
#'
#' # Create `strollur::strollur` object
#' data <- strollur::new_dataset("example")
#'
#' # Add FASTQ data `strollur::strollur` object
#' strollur::add(data, table, type = "fastq")
#'
#' # Write `strollur::strollur` objects sequence quality data to file
#' strollur::write_quality(data, tempfile())
#'
#' @return name of sequence quality file
#' @export
write_quality <- function(data, filename = NULL) {
  if (!inherits(data, "strollur")) {
    stop("data must be a strollur object.")
  }

  if (is.null(filename)) {
    filename <- names(data, "dataset")
    if (filename == "") {
      .abort_no_name()
    }
    filename <- paste0(filename, ".qual")
  }

  quality <- xdev_report(data, type = "quality")

  # data contains sequences
  if (nrow(quality) != 0) {
    clean_scores <- sapply(
      quality$quality_score,
      function(x) paste(x, collapse = ", ")
    )

    # Format the data
    formatted_lines <- paste0(
      ">", quality$sequence_name,
      "\n", clean_scores
    )

    # Write to a text file
    writeLines(formatted_lines, filename)
    return(filename)
  }

  "no_quality_data"
}
