#' @title export
#' @description
#' The export function will create a list containing the data in the dataset.
#' @param dataset a 'dataset' object
#' @param tags a vector of strings containing the items you wish to export.
#' Options are 'sequence_data', 'bin_data', 'metadata',
#' 'references', 'sequence_tree', 'sample_tree', 'alignment_report',
#' 'contigs_assembly_report' and 'chimera_report'. By default, everything is
#'  exported.
#'
#' @examples
#'
#' miseq <- miseq_sop_example()
#' table <- export(miseq)
#'
#' @return A list of data.frames
#' @seealso [dataset$export()]
#' @export
export <- function(dataset, tags = NULL) {
  dataset$export(tags)
}
