#' @title export
#' @description
#' The export function will create a list containing the data in the dataset.
#' @param dataset a 'sequence_data' object
#' @param tags a vector of strings containing the items you wish to export.
#' Options are 'sequence_data', 'bin_data', 'metadata',
#' 'references', 'sequence_tree', 'sample_tree', 'alignment_report',
#' 'contigs_assembly_report'. By default, everything is exported.
#'
#' @examples
#'
#' miseq <- miseq_sop_example()
#' table <- export(miseq)
#'
#' @return A list of data.frames
#' @seealso [sequence_data$export()]
#' @export
export <- function(dataset, tags = NULL) {
  dataset$export(tags)
}
