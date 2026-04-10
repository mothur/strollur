#' @title write_mothur_design
#' @description
#' Write a mothur formatted
#' \href{https://mothur.org/wiki/design_file/}{design file}
#'
#' @param data A `strollur` object
#' @param filename a string containing the name of the output file. Default =
#' 'dataset_name'.design
#'
#' @examples
#'
#' miseq <- miseq_sop_example()
#' write_mothur_design(miseq, tempfile())
#'
#' @return name of design file
#' @export
write_mothur_design <- function(data, filename = NULL) {
  if (!inherits(data, "strollur")) {
    stop("data must be a strollur object.")
  }

  if (is.null(filename)) {
    filename <- names(data, "dataset")
    if (filename == "") {
      .abort_no_name()
    }
    filename <- paste0(filename, ".design")
  }

  df <- report(data, "sample_assignments")

  if (nrow(df) != 0) {
    # write table
    readr::write_tsv(df, filename, col_names = TRUE)
    return(filename)
  }

  "no_design_data"
}
