#' @title load
#' @description
#' The load function will create a 'sequence_data' object from a file.
#' @param file a string containing the file name.
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
  dataset$data <- new(Dataset, "", 1)
  dataset$deserialize()
}
