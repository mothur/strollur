#' @title is_equal
#' @description
#' Determine if two
#' \href{https://mothur.org/strollur/reference/strollur.html}{strollur object}s
#' are equal.
#'
#' @param data, a
#'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur object}
#' @param data2, a
#'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur object}
#'
#' @examples
#'
#' miseq <- strollur::miseq_sop_example()
#'
#' data <- strollur::copy_dataset(miseq)
#'
#' strollur::is_equal(miseq, data)
#'
#' @returns a logical
#' @export
is_equal <- function(data, data2) {
  if (!inherits(data, "strollur")) {
    stop("data must be a strollur object.")
  }

  if (!inherits(data2, "strollur")) {
    stop("data2 must be a strollur object.")
  }

  data$is_equal(data2)
}
