#' @title read_mothur_list
#' @description
#' Read a mothur formatted \href{https://mothur.org/wiki/list_file/}{list file}
#' @param list file name (required)
#' @examples
#'
#' # You can add your otu assignments to the your data set using the following:
#'
#' otu_data <- read_mothur_list(rdataset_example("final.opti_mcc.list"))
#'
#' dataset <- sequence_data$new()
#' dataset$assign_bins(otu_data)
#'
#' @return A data.frame containing the sequence otu assignments
#' @export
read_mothur_list <- function(list) {
  if (!file.exists(list)) {
    abort_nonexistant_file(list)
  }

  df <- readr::read_table(
    file = list, col_names = TRUE,
    show_col_types = FALSE
  )

  # remove label and numOtus columns
  df <- df[, -c(1, 2)]
  otu_names <- names(df)

  df <- apply(df, 2, split_at_char)

  otu_assignments <- c()
  sequence_names <- c()

  i <- 1
  for (otu in df) {
    sequence_names <- c(sequence_names, otu)
    otu_assignments <- c(
      otu_assignments,
      rep(otu_names[i], length(otu))
    )
    i <- i + 1
  }

  data.frame(
    bin_names = otu_assignments,
    sequence_names = sequence_names
  )
}
