#' @title write_mothur_list
#' @description
#' Write mothur formatted \href{https://mothur.org/wiki/list_file/}{list files}
#'
#' @param data A 'dataset' object
#' @param file_root a string containing the root name of the output file.
#' Default = 'dataset_name'. Resulting in output files
#' 'dataset_name'.bin_type'.list.
#'
#' @examples
#'
#' miseq <- miseq_sop_example()
#' write_mothur_list(miseq)
#'
#' @return vector containing the names of the files created
#' @export
write_mothur_list <- function(data, file_root = NULL) {
  # check type
  if (class(data)[1] != "dataset") {
    .abort_incorrect_type("dataset", data)
  }

  if (is.null(file_root)) {
    file_root <- name(data, "dataset")
    if (file_root == "") {
      .abort_no_name()
    }
  }

  bin_types <- data$get_bin_types()
  outputs <- c()

  for (type in bin_types) {
    df <- get_list(data, type)

    if (nrow(df) != 0) {
      output_file <- paste0(file_root, ".", type, ".list")
      outputs <- c(outputs, output_file)

      header <- paste("label",
        paste0("num", type, "s"),
        paste(unique(df[[1]]), collapse = "\t"),
        collapse = "\t"
      )

      t <- tapply(df[[2]], df[[1]], paste, collapse = ",")

      bins <- paste("1",
        length(t),
        paste(t, collapse = "\t"),
        collapse = "\t"
      )

      lines <- c(header, bins)
      readr::write_lines(lines, output_file)
    }
  }
  outputs
}
