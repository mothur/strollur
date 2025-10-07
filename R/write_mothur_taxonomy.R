#' @title write_mothur_taxonomy
#' @description
#' Write a mothur formatted
#' \href{https://mothur.org/wiki/taxonomy_file/}{taxonomy file}
#'
#' @param data A 'dataset' object
#' @param filename a string containing the name of the output file. Default =
#' 'dataset_name'.taxonomy
#'
#' @examples
#'
#' miseq <- miseq_sop_example()
#' write_mothur_taxonomy(miseq)
#'
#' @return name of taxonomy file
#' @export
write_mothur_taxonomy <- function(data, filename = NULL) {
  # check type
  if (class(data)[1] != "dataset") {
    abort_incorrect_type("dataset", data)
  }

  if (is.null(filename)) {
    filename <- data$get_dataset_name()
    if (filename == "") {
      abort_no_name()
    }
    filename <- paste0(filename, ".taxonomy")
  }

  df <- data$get_sequence_taxonomy_report()

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

  return("no_sequence_taxonomy")
}
