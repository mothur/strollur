#' @title read_mothur_cons_taxonomy
#' @description
#' Read a mothur formatted
#' \href{https://mothur.org/wiki/constaxonomy_file/}{constaxonomy file}
#' @param taxonomy file name (required)
#' @examples
#'
#' # You can add the otu assignments and bin taxonomies to the your data set
#' # using the following:
#'
#' otu_data <- read_mothur_cons_taxonomy(rdataset_example(
#'   "final.cons.taxonomy"
#' ))
#'
#' dataset <- sequence_data$new()
#' dataset$assign_bins(otu_data)
#' dataset$assign_bin_taxonomy(otu_data)
#'
#' @return A data.frame containing the bin names, bin abundances and bin
#' taxonomies.
#' @export
read_mothur_cons_taxonomy <- function(taxonomy) {
  if (!file.exists(taxonomy)) {
    abort_nonexistant_file(taxonomy)
  }

  df <- readr::read_table(
    file = taxonomy, col_names = TRUE,
    show_col_types = FALSE
  )

  names(df) <- c("bin_names", "abundances", "taxonomies")
  df
}
