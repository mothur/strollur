#' @title write_mothur_shared
#' @description
#' Write mothur formatted
#' \href{https://mothur.org/wiki/shared_file/}{shared files}
#'
#' @param data A `strollur` object
#' @param file_root a string containing the root name of the output file.
#' Default = 'dataset_name'. Resulting in output files
#' 'dataset_name'.bin_type'.shared.
#'
#' @examples
#'
#' miseq <- miseq_sop_example()
#' write_mothur_shared(miseq, tempfile())
#'
#' @return vector containing the names of the files created
#' @export
write_mothur_shared <- function(data, file_root = NULL) {
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
    df <- abundance(
      data = data, type = "bins",
      bin_type = type, by_sample = TRUE
    )

    if (nrow(df) != 0) {
      num_cols <- ncol(df)

      # you need bin_names, and abundances by sample for a shared file
      if (num_cols > 2) {
        output_file <- paste0(file_root, ".", type, ".shared")
        outputs <- c(outputs, output_file)

        # remove treatment columns if present
        if (num_cols == 4) {
          df <- df[, -c(4)]
        }

        col_names <- names(df)

        # create lines for samples
        df <- df |>
          pivot_wider(
            names_from = col_names[1],
            values_from = col_names[2],
            id_cols = col_names[3],
            values_fill = 0
          )

        number_of_bins <- ncol(df) - 1
        num_samples <- nrow(df)

        # add label column before samples
        # To fix build warnings, we have to make bindings for:
        # label, num_bins, samples
        label <- num_bins <- samples <- NULL

        df <- df |>
          mutate(label = rep(1, num_samples)) |>
          relocate(label, .before = samples)

        # add numotus column after samples
        df <- df |>
          mutate(num_bins = rep(number_of_bins, num_samples)) |>
          relocate(num_bins, .after = samples)

        readr::write_tsv(df, output_file)
      }
    }
  }
  outputs
}
