#' @title load
#' @description
#' The load function will create a 'sequence_data' object from an RDS file.
#' @param file a string containing the .rds file name.
#' @examples
#'
#' dataset <- load(rdataset_example("miseq_sop.rds"))
#' dataset
#'
#' @return A 'sequence_data' object
#' @export
load <- function(file) {
  if (!file.exists(file)) {
    abort_nonexistant_file(file)
  }

  dataset <- readRDS(file)
  dataset$data <- new_dataset("", 1)
  load_dataset(dataset$data, dataset$raw)
  dataset$raw <- NULL
  dataset
}
