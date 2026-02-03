#' @title save_dataset
#' @description
#' The save_dataset function will save the \link{dataset} object to file.
#'
#' @param dataset a \link{dataset} object
#'
#' @param file a string containing the file name.
#' @examples
#'
#' dataset <- read_mothur(
#'   fasta = rdataset_example("final.fasta.gz"),
#'   count = rdataset_example("final.count_table.gz"),
#'   taxonomy = rdataset_example("final.taxonomy.gz"),
#'   design = rdataset_example("mouse.time.design"),
#'   otu_list = rdataset_example("final.opti_mcc.list.gz"),
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
