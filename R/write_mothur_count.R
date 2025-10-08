#' @title write_mothur_count
#' @description
#' Write a mothur formatted
#' \href{https://mothur.org/wiki/count_file/}{count file}
#'
#' @param data A 'dataset' object
#' @param filename a string containing the name of the output file. Default =
#' 'dataset_name'.count_table
#'
#' @examples
#'
#' miseq <- miseq_sop_example()
#' write_mothur_count(miseq)
#'
#' @return name of count file
#' @export
write_mothur_count <- function(data, filename = NULL) {
  # check type
  if (class(data)[1] != "dataset") {
    abort_incorrect_type("dataset", data)
  }

  if (is.null(filename)) {
    filename <- data$get_dataset_name()
    if (filename == "") {
      abort_no_name()
    }
    filename <- paste0(filename, ".count_table")
  }

  df <- data$get_count_table()

  if (nrow(df) != 0) {
    # treatment assignments, remove treatment column
    if (ncol(df) == 4) {
      df <- df[, -c(4)]
    }

    # has samples
    if (ncol(df) == 3) {
      samples <- sort(unique(df[[3]]))
      compressed_header <- paste0("#Compressed Format: groupIndex,abunda",
        "nce. For example 1,6 would mean the ",
        "read has an abundance of 6 for group ",
        samples[1], ".\n#",
        collapse = ""
      )

      # print compressed header name assignments ie. 1,F0D3\t
      compressed_header <- paste0(compressed_header,
        paste(seq_along(samples), samples,
          sep = ",", collapse = "\t"
        ),
        "\nRepresentative_Sequence\ttotal\t",
        paste(samples, collapse = "\t"),
        collapse = ""
      )
      # write headers
      writeLines(compressed_header, filename)

      # Create a mapping from unique sample names to a new index
      # The factor levels are used to create the index
      sample_map <- data.frame(
        original_sample = samples,
        sample_index = as.integer(as.factor(samples))
      )

      # Join the original dataframe with the new mapping to replace
      # sample names
      df <- left_join(df, sample_map,
        by = c("samples" = "original_sample")
      )

      # Group by 'ids' and summarize the 'sample_index' and 'abunds'
      df <- df %>%
        group_by(sequence_names) %>%
        summarise(
          total = sum(abundances),
          sample_abund = paste0(sample_index, ",", abundances,
            collapse = " "
          )
        )

      # write table
      readr::write_tsv(df, filename, append = TRUE, col_names = FALSE)
    } else {
      names(df) <- c("Representative_Sequence", "total")

      # write table data
      readr::write_tsv(df, filename)
    }

    return(filename)
  }

  return("no_sequence_data")
}
