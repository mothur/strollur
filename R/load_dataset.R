#' @title Load strollur object from .rds file
#' @description The load_dataset function will create a
#'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#'   from an RDS file.
#'
#' @param file a string containing the .rds file name.
#' @examples
#'
#' data <- load_dataset(strollur_example("miseq_sop.rds"))
#' data
#'
#' @seealso [save_dataset()]
#' @return a
#'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#' @export
load_dataset <- function(file) {
  if (!file.exists(file)) {
    .abort_nonexistant_file(file)
  }

  # R6 strollur class. back end c++ Dataset class needs to be built
  # using xint_deserialize_dobject.
  dataset <- readRDS(file)

  if (!inherits(dataset, "strollur")) {
    stop("The .rds file must contain a strollur object.")
  }

  dataset_version <- dataset$get_version()
  current_version <- new_dataset()$get_version()

  # check attributes for valid version
  if (current_version != dataset_version) {
    if (.import_compatible(dataset_version)) {
      message <- paste0(
        "The table was created with strollur version ", dataset_version,
        ", which is compatible with the current strollur version: ",
        current_version, ". Converting and importing."
      )
      cli::cli_alert(message)
    } else {
      message <- paste0(
        "Unable to create 'strollur' object. ",
        "The table was created with strollur version ", dataset_version,
        ", which is not compatible with the current strollur version: ",
        current_version, "."
      )
      cli::cli_abort(message)
    }
  }

  xint_deserialize_dobject(dataset)
  dataset
}
