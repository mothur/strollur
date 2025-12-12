#' @title write_mothur_cons_taxonomy
#' @description
#' Write a mothur formatted
#' \href{https://mothur.org/wiki/constaxonomy_file/}{cons_taxonomy file}
#'
#' @param data A 'dataset' object
#' @param file_root a string containing the root name of the output file.
#' Default = 'dataset_name'. Resulting in output files
#' 'dataset_name'.bin_type'.cons.taxonomy.
#'
#' @examples
#'
#' miseq <- miseq_sop_example()
#' write_mothur_cons_taxonomy(miseq)
#'
#' @return vector containing the names of the files created
#' @export
write_mothur_cons_taxonomy <- function(data, file_root = NULL) {
  # check type
  if (class(data)[1] != "dataset") {
    .abort_incorrect_type("dataset", data)
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
    df <- report(data, "bin_taxonomy", type)

    if (nrow(df) != 0) {
      # confidence scores are available
      if (ncol(df) == 4) {
        # create string like "Bacteria(100);"
        df[[3]] <- paste0(df[[3]], "(", df[[4]], ");")
      }

      df <- data.frame(
        bin_name = unique(df[[1]]),
        bin_abundance = abundance(
          data = data,
          type = "bins", bin_type = type
        )[[2]],
        bin_taxonomy = tapply(df[[3]], df[[1]],
          paste,
          collapse = ""
        )
      )

      output_file <- paste0(file_root, ".", type, ".cons.taxonomy")
      outputs <- c(outputs, output_file)

      readr::write_tsv(df, output_file, escape = "none")
    }
  }
  outputs
}
