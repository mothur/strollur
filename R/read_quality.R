#' @title Read a file containing sequence quality scores
#'
#' @description
#' Read a file containing sequence quality scores
#'
#' @param quality string, containing the quality score file name.
#' @param path string, optional path to quality score file.
#'
#' @examples
#'
#' quality_data <- strollur::read_quality(strollur_example("tiny.qual"))
#'
#' # quality_data is a data.frame.
#' # To access the names of the sequences in the file, run the following:
#'
#' quality_data$sequence_name
#'
#' # To access the quality scores in the file, run the following:
#'
#' quality_data$quality_score
#'
#' @return A data.frame containing the quality score sequence data
#' @export
read_quality <- function(quality, path = NULL) {
  # if path is set, adjust file name
  if (!is.null(path)) {
    quality <- file.path(path, basename(quality))
  }

  if (!file.exists(quality)) {
    .abort_nonexistant_file(quality)
  }

  lines <- readLines(quality)
  is_sequence_name <- grepl("^>", lines)

  # remove '>'
  sequence_names <- sub("^>", "", lines[is_sequence_name])
  # remove comments
  sequence_names <- gsub(" .*$", "", sequence_names)

  scores_raw <- lines[!is_sequence_name]
  scores_list <- lapply(scores_raw, function(x) {
    # Remove trailing/leading spaces and split by comma
    as.numeric(strsplit(trimws(x), ",\\s*")[[1]])
  })

  data.frame(
    sequence_name = sequence_names,
    quality_score = I(scores_list)
  )
}
