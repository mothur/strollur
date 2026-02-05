#' @title read_mothur_rabund
#' @description
#' Read a mothur formatted
#' \href{https://mothur.org/wiki/rabund_file/}{rabund file}
#' @param rabund file name (required)
#' @examples
#'
#' # You can add your otu assignments to the your data set using the following:
#'
#' # read rabund file into data.frame
#' otu_data <- read_mothur_rabund(
#'   rabund =
#'     strollur_example("final.opti_mcc.rabund")
#' )
#'
#  #create a new empty dataset
#' data <- new_dataset()
#'
#' # assign abundance only 'otu' bins
#' assign(data = data, table = otu_data, type = "bins", bin_type = "otu")
#'
#' @return A data.frame containing the sequence otu assignments
#' @export
read_mothur_rabund <- function(rabund) {
  if (!file.exists(rabund)) {
    .abort_nonexistant_file(rabund)
  }

  df <- readr::read_table(
    file = rabund, col_names = TRUE,
    show_col_types = FALSE
  )

  df <- df[, -c(1, 2)]

  data.frame(
    bin_names = names(df),
    abundances = t(df)
  )
}
