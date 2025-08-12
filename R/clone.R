#' @title clone
#' @description
#' The clone function will create a deep copy of 'sequence_data' object.
#' @param dataset a 'sequence_data' object
#' @examples
#' # For dataset's including sequence data:
#'
#' dataset <- read_mothur(
#'   fasta = rdataset_example("final.fasta"),
#'   count = rdataset_example("final.count_table"),
#'   taxonomy = rdataset_example("final.taxonomy"),
#'   design = rdataset_example("mouse.time.design"),
#'   list = rdataset_example("final.opti_mcc.list"),
#'   dataset_name = "miseq_sop"
#' )
#'
#' copy_dataset <- clone(dataset)
#' copy_dataset
#'
#' @return A 'sequence_data' object
#' @export
clone <- function(dataset) {
  sequence_data$new(dataset = dataset)
}
