#' @title save
#' @description
#' The save function will save the 'sequence_data' object to file.
#' @param dataset a 'sequence_data' object
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
#' save(dataset, "miseq_sop.rds")
#'
#' @return A file containing the 'sequence_data' object
#' @export
save <- function(dataset, file) {
  if (class(dataset)[1] != "sequence_data") {
    abort_incorrect_type("sequence_data", class(dataset)[1])
  }

  dataset$serialize()
  saveRDS(dataset, file = file)
}
