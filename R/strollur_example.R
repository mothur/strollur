#' Get path to strollur example
#'
#' @title strollur_example
#' @description
#' strollur comes bundled with some example files in its `inst/extdata`
#' directory. This function make them easy to access them.
#'
#' @param file Name of file.
#' @examples
#' strollur_example()
#' strollur_example("final.fasta.gz")
#'
#' @return string, Full path to example files
#' @export
strollur_example <- function(file = NULL) {
  path <- ""

  if (is.null(file)) {
    path <- system.file("extdata", package = "strollur")
  } else {
    path <- system.file("extdata", file,
      package = "strollur",
      mustWork = TRUE
    )
  }
  path
}
