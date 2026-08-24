#' @title export_dataset
#' @description
#' Export all data from a
#' \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object.
#'
#' @param data a
#'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#'
#' @examples
#'
#' dataset <- strollur::new_dataset("my_dataset")
#' strollur::export_dataset(dataset)
#'
#' @return list, containing the data in the strollur::strollur object
#' @export
export_dataset <- function(data) {
  if (!inherits(data, "strollur")) {
    stop("data must be a strollur object.")
  }

  results <- xdev_export_dataset(data)

  if (!is.null(data$sequence_tree)) {
    results[["sequence_tree"]] <- data$sequence_tree
  }

  if (!is.null(data$sample_tree)) {
    results[["sample_tree"]] <- data$sample_tree
  }

  attr(
    results,
    "strollur_version"
  ) <- data[[".__enclos_env__"]]$private$version
  attr(results, "dataset_name") <- names(data, type = "dataset")

  results
}
