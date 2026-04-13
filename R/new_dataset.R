#' @title new_dataset
#' @description
#' Create a new \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#'
#' @param dataset_name string, a string containing the dataset name.
#' Default = ""
#' @param processors integer, number of cores to use during summary functions.
#' Default = all available
#' @examples
#'
#' data <- new_dataset()
#'
#' # to create a new dataset named "soil" and allow for all available
#' # processors during summary functions, run the following:
#'
#' data <- new_dataset(dataset_name = "soil")
#'
#' # to create a new dataset named "soil" and allow for 2
#' # processors during summary functions, run the following:
#'
#' data <- new_dataset(dataset_name = "soil", processors = 2)
#'
#' @returns a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#' @export
#' @seealso The 'new' method in the \href{https://mothur.org/strollur/reference/strollur.html}{strollur} class
new_dataset <- function(dataset_name = "", processors = NULL) {
  if(is.null(processors)) {
    processors <- parallelly::availableCores()
  }
  strollur$new(dataset_name, processors)
}