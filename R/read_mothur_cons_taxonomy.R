#' @title read_mothur_cons_taxonomy
#' @description
#' Read a mothur formatted
#' \href{https://mothur.org/wiki/constaxonomy_file/}{cons_taxonomy file}
#' @param taxonomy file name, a mothur
#' \href{https://mothur.org/wiki/constaxonomy_file/}{consensus taxonomy file}.
#'  The cons_taxonomy file is created by
#' \href{https://mothur.org/wiki/classify.otu/}{classify.otu}.
#'
#' @examples
#'
#' # You can add the otu assignments and bin taxonomies to the your data set
#' # using the following:
#'
#' # read mothur's consensus taxonomy file into a data.frame
#' otu_data <- read_mothur_cons_taxonomy(rdataset_example(
#'   "final.cons.taxonomy"
#' ))
#'
#' data <- new_dataset()
#'
#' # assign abundance only 'otu' bins
#' assign(data = data, table = otu_data, type = "bins", bin_type = "otu")
#'
#' # assign consensus taxonomies to 'otu' bins
#' assign(
#'   data = data, table = otu_data,
#'   type = "bin_taxonomy", bin_type = "otu"
#' )
#'
#' @return A data.frame containing the bin names, bin abundances and bin
#' taxonomies.
#' @export
read_mothur_cons_taxonomy <- function(taxonomy) {
  if (!file.exists(taxonomy)) {
    .abort_nonexistant_file(taxonomy)
  }

  df <- readr::read_table(
    file = taxonomy, col_names = TRUE,
    show_col_types = FALSE
  )

  names(df) <- c("bin_names", "abundances", "taxonomies")
  df
}
