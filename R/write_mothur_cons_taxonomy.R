#' @title write_mothur_cons_taxonomy
#' @description
#' Write a mothur formatted
#' \href{https://mothur.org/wiki/constaxonomy_file/}{cons_taxonomy file}
#'
#' @param data A 'dataset' object
#' @param filename a string containing the name of the output file. Default =
#' 'dataset_name'.'bin_type'.cons.taxonomy
#'
#' @examples
#'
#' miseq <- miseq_sop_example()
#' write_mothur_cons_taxonomy(miseq)
#'
#' @return vector containing the names of the files created
#' @export
write_mothur_cons_taxonomy <- function(data, filename = NULL) {
  # check type
  if (class(data)[1] != "dataset") {
    abort_incorrect_type("dataset", data)
  }

  if (is.null(filename)) {
    filename <- data$get_dataset_name()
    if (filename == "") {
      abort_no_name()
    }
  }

  bin_types <- data$get_bin_types()
  outputs <- c()

  for (type in bin_types) {
    df <- data$get_bin_taxonomy_report(type)

    if (nrow(df) != 0) {
      # confidence scores are available
      if (ncol(df) == 4) {
        # create string like "Bacteria(100);"
        df[[3]] <- paste0(df[[3]], "(", df[[4]], ");")
      }

      df <- data.frame(
        bin_name = unique(df[[1]]),
        bin_abundance = get_rabund_vector(data$data, type),
        bin_taxonomy = tapply(df[[3]], df[[1]],
          paste,
          collapse = ""
        )
      )

      output_file <- paste0(filename, ".", type, ".cons.taxonomy")
      outputs <- c(outputs, output_file)

      readr::write_tsv(df, output_file, escape = "none")
    }
  }
  outputs
}
