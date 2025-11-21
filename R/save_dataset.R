#' @title save_dataset
#' @description
#' The save_dataset function will save the \link{dataset} object to file.
#' @param dataset a \link{dataset} object
#' @param file a string containing the file name.
#' @examples
#'
#' dataset <- read_mothur(
#'   fasta = rdataset_example("final.fasta"),
#'   count = rdataset_example("final.count_table"),
#'   taxonomy = rdataset_example("final.taxonomy"),
#'   design = rdataset_example("mouse.time.design"),
#'   otu_list = rdataset_example("final.opti_mcc.list"),
#'   dataset_name = "miseq_sop"
#' )
#'
#' save_dataset(dataset, "miseq_sop.rds")
#'
#' @return A file containing the 'dataset' object
#' @export
save_dataset <- function(dataset, file) {
  if (class(dataset)[1] != "dataset") {
    .abort_incorrect_type("dataset", dataset)
  }

  xint_serialize_dobject(dataset)
  saveRDS(dataset, file = file)
  file
}
