#' @title load_dataset
#' @description
#' The load_dataset function will create a \link{dataset} object from an RDS
#' file.
#' @param file a string containing the .rds file name.
#' @examples
#'
#' data <- load_dataset(strollur_example("miseq_sop.rds"))
#' data
#'
#' @return a \link{dataset} object
#' @export
load_dataset <- function(file) {
  if (!file.exists(file)) {
    .abort_nonexistant_file(file)
  }

  dataset <- readRDS(file)
  xint_deserialize_dobject(dataset)
  dataset
}
