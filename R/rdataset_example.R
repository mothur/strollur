#' Get path to rdataset example
#'
#' rdataset comes bundled with some example files in its `inst/extdata`
#' directory. This function make them easy to access.
#'
#' @param file Name of file.
#' @export
#' @examples
#' rdataset_example()
#' rdataset_example("final.fasta")
rdataset_example <- function(file = NULL) {
  path <- ""

  if (is.null(file)) {
    path <- system.file("extdata", package = "rdataset")
  } else {
    path <- system.file("extdata", file,
      package = "rdataset",
      mustWork = TRUE
    )
  }
  return(path)
}
