#' @title write_mothur_rabund
#' @description
#' Write mothur formatted
#' \href{https://mothur.org/wiki/rabund_file/}{rabund files}
#'
#' @param data A `strollur` object
#' @param file_root a string containing the root name of the output file.
#' Default = 'dataset_name'. Resulting in output files
#' 'dataset_name'.bin_type'.rabund.
#'
#' @examples
#'
#' miseq <- miseq_sop_example()
#' write_mothur_rabund(miseq, tempfile())
#'
#' @return vector containing the names of the files created
#' @export
write_mothur_rabund <- function(data, file_root = NULL) {
  # check type
  if (class(data)[1] != "strollur") {
    .abort_incorrect_type("strollur", data)
  }

  if (is.null(file_root)) {
    file_root <- names(data, "dataset")
    if (file_root == "") {
      .abort_no_name()
    }
  }

  bin_types <- data$get_bin_types()
  outputs <- c()

  for (type in bin_types) {
    df <- abundance(data = data, type = "bins", bin_type = type)

    if (nrow(df) != 0) {
      output_file <- paste0(file_root, ".", type, ".rabund")
      outputs <- c(outputs, output_file)

      header <- paste("label",
        paste0("num", type, "s"),
        paste(unique(df[[1]]), collapse = "\t"),
        collapse = "\t"
      )

      bins <- paste("1",
        nrow(df),
        paste(df[[2]], collapse = "\t"),
        collapse = "\t"
      )

      lines <- c(header, bins)
      readr::write_lines(lines, output_file)
    }
  }
  outputs
}
