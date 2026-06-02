#' @title new_dataset
#' @description Create a new
#' \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#'
#' @param dataset_name string, a string containing the dataset name.
#' Default = ""
#' @examples
#'
#' data <- new_dataset()
#'
#' # to create a new dataset named "soil", run the following:
#'
#' data <- new_dataset(dataset_name = "soil")
#'
#' @returns a
#'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#' @export
#' @seealso The 'new' method in the
#'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur} class
new_dataset <- function(dataset_name = "") {
  strollur$new(dataset_name)
}
