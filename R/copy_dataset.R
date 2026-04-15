#' @title copy_dataset
#' @description
#' Create a new
#' \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#' from an existing dataset.
#'
#' @param data, a
#'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#' @examples
#'
#' miseq <- miseq_sop_example()
#'
#' # to create a new dataset that is a copy of miseq
#'
#' data <- copy_dataset(miseq)
#'
#' @returns a
#'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#' @seealso The 'new' method in the
#'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur} class
#' @export
copy_dataset <- function(data) {
  if (!inherits(data, "strollur")) {
    stop("data must be a strollur object.")
  }

  copy <- strollur$new(
    xdev_names(data,
      type = "dataset"
    ),
    parallelly::availableCores(),
    data
  )
}
