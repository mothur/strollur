#' @title clone
#' @description
#' The clone function will create a deep copy of 'dataset' object. Note:
#' 'dataset' is an R6 class with reference semantics and requires a deep
#' copy.
#' @param data a 'dataset' object
#' @examples
#' # For dataset's including sequence data:
#'
#' data <- read_mothur(
#'   fasta = rdataset_example("final.fasta"),
#'   count = rdataset_example("final.count_table"),
#'   taxonomy = rdataset_example("final.taxonomy"),
#'   design = rdataset_example("mouse.time.design"),
#'   otu_list = rdataset_example("final.opti_mcc.list"),
#'   dataset_name = "miseq_sop"
#' )
#'
#' copy_dataset <- clone(data)
#' copy_dataset
#'
#' @return A 'dataset' object
#' @export
clone <- function(data) {
  dataset$new(dataset = data)
}
