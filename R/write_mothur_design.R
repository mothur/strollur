#' @title write_mothur_design
#' @description
#' Write a mothur formatted
#' \href{https://mothur.org/wiki/design_file/}{design file}
#'
#' @param data A 'dataset' object
#' @param filename a string containing the name of the output file. Default =
#' 'dataset_name'.design
#'
#' @examples
#'
#' miseq <- miseq_sop_example()
#' write_mothur_design(miseq)
#'
#' @return name of design file
#' @export
write_mothur_design <- function(data, filename = NULL) {
  # check type
  if (class(data)[1] != "dataset") {
    abort_incorrect_type("dataset", data)
  }

  if (is.null(filename)) {
    filename <- data$get_dataset_name()
    if (filename == "") {
      abort_no_name()
    }
    filename <- paste0(filename, ".design")
  }

  df <- data$get_sample_treatment_assignments()

  if (nrow(df) != 0) {
    # write table
    readr::write_tsv(df, filename, col_names = TRUE)
    return(filename)
  }

  return("no_design_data")
}
