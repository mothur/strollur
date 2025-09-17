#' @title read_mothur_taxonomy
#' @description
#' Read a mothur formatted
#' \href{https://mothur.org/wiki/taxonomy_file/}{taxonomy file}
#' @param taxonomy file name (required)
#' @examples
#'
#' # You can add the sequences and their taxonomies to the your data set
#' # using the following:
#'
#' classification_data <- read_mothur_taxonomy(rdataset_example(
#'   "final.taxonomy"
#' ))
#'
#' dataset <- sequence_data$new()
#' dataset$assign_sequence_taxonomy(classification_data)
#'
#' @return A data.frame containing the sequences names and sequences taxonomies.
#' @export
read_mothur_taxonomy <- function(taxonomy) {
  if (!file.exists(taxonomy)) {
    abort_nonexistant_file(taxonomy)
  }

  df <- readr::read_table(
    file = taxonomy, col_names = FALSE,
    show_col_types = FALSE
  )

  names(df) <- c("sequence_names", "taxonomies")
  df
}
