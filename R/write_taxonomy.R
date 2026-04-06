#' @title write_taxonomy
#' @description
#' Write a 2 column
#' \href{https://mothur.org/wiki/taxonomy_file/}{taxonomy file}
#'
#' @param data A `strollur` object
#' @param filename a string containing the name of the output file. Default =
#' 'dataset_name'.taxonomy
#'
#' @examples
#'
#' miseq <- miseq_sop_example()
#' write_taxonomy(miseq, tempfile())
#'
#' @return name of taxonomy file
#' @export
write_taxonomy <- function(data, filename = NULL) {
  # check type
  if (class(data)[1] != "strollur") {
    .abort_incorrect_type("strollur", data)
  }

  if (is.null(filename)) {
    filename <- names(data, "dataset")
    if (filename == "") {
      .abort_no_name()
    }
    filename <- paste0(filename, ".taxonomy")
  }

  df <- report(data, "sequence_taxonomy")

  if (nrow(df) != 0) {
    # confidence scores are available
    if (ncol(df) == 4) {
      # create string like "Bacteria(100);"
      df[[3]] <- paste0(df[[3]], "(", df[[4]], ");")
    }

    df <- data.frame(
      unique(df[[1]]),
      tapply(df[[3]], df[[1]], paste, collapse = "")
    )

    readr::write_tsv(df, filename, escape = "none", col_names = FALSE)
    return(filename)
  }

  "no_sequence_taxonomy"
}
