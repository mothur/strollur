#' @title read_mothur_rabund
#' @description
#' Read a mothur formatted
#' \href{https://mothur.org/wiki/rabund_file/}{rabund file}
#' @param rabund file name (required)
#' @examples
#'
#' # You can add your otu assignments to the your data set using the following:
#'
#' otu_data <- read_mothur_rabund(rdataset_example("final.opti_mcc.rabund"))
#'
#' dataset <- sequence_data$new()
#' dataset$assign_bins(otu_data$bin_id, otu_data$abundance)
#'
#' @return A data.table containing the sequence otu assignments
#' @export
read_mothur_rabund <- function(rabund) {
  if (!file.exists(rabund)) {
    abort_nonexistant_file(rabund)
  }

  df <- readr::read_table(
    file = rabund, col_names = FALSE,
    show_col_types = FALSE
  )
  num_otus <- df[[2]]
  df <- df[, -c(1, 2)]

  otu_assignments <- paste0("otu", c(1:num_otus))
  abundance <- transpose(df)[[1]]

  data.table(
    bin_id = otu_assignments,
    abundance = abundance
  )
}
