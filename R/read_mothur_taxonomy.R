#' @title read_mothur_taxonomy
#' @description
#' Read a mothur formatted
#' \href{https://mothur.org/wiki/taxonomy_file/}{taxonomy file}
#' @param taxonomy file name. a mothur
#' \href{https://mothur.org/wiki/taxonomy_file/}{taxonomy file}, created by
#' \href{https://mothur.org/wiki/classify.seqs/}{classify.seqs}
#' @examples
#'
#' # You can add the sequences and their taxonomies to the your data set
#' # using the following:
#'
#' # read mothur's taxonomy file into a data.frame
#' classification_data <- read_mothur_taxonomy(strollur_example(
#'   "final.taxonomy.gz"
#' ))
#'
#' # create a new empty `strollur` object
#' data <- new_dataset()
#'
#' # assign sequence classifications
#' assign(data = data, table = classification_data, type = "sequence_taxonomy")
#'
#' @return A data.frame containing the sequences names and sequences taxonomies.
#' @export
read_mothur_taxonomy <- function(taxonomy) {
  if (!file.exists(taxonomy)) {
    .abort_nonexistant_file(taxonomy)
  }

  df <- readr::read_table(
    file = taxonomy, col_names = FALSE,
    show_col_types = FALSE
  )

  names(df) <- c("sequence_names", "taxonomies")
  df
}
