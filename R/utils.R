# =========================================================================== #
#' @title remove_file
#' @description Remove file, if it exists
#' @param filename, String containing name of file to remove
remove_file <- function(filename) {
  if (file.exists(filename)) {
    file.remove(filename)
  }
}
# =========================================================================== #
#' @title sort_dataframe
#' @description Sort dataframe
#' @param data, the data.frame to be sorted
#' @param order, vector containing the order desired
#' @param named_col, name of column in data.frame to match order
#'
#' # sort results alphabetically
#'
#' miseq <- miseq_sop_example()
#'
#' sequence_names <- names(miseq)
#'
#' fasta <- report(miseq, type = fasta)
#'
#' sorted_fasta <- sort_dataframe(fasta,
#'                                order = sort(sequence_names),
#'                                named_col = "sequence_names")
#'
#' @return sorted data.frame
#' @export
sort_dataframe <- function(data, order, named_col) {
  sorted <- data[order(match(data[[named_col]], order)), ]
}
# =========================================================================== #
#' @title .split_white_space
#' @description Split string at white space
#' @param line, String containing data to split
#'
#' @returns A vector of Strings
#' @keywords internal
#' @noRd
.split_white_space <- function(line) {
  words <- strsplit(line, "\\s")[[1]]
  words <- words[nzchar(x = words)]
}
# =========================================================================== #
#' @title .split_at_char
#' @description Split string by deliminating character
#' @param line, String containing data to split
#' @param delim, Character to split data by
#' @returns A vector of Strings
#' @keywords internal
#' @noRd
.split_at_char <- function(line, delim = ",") {
  words <- strsplit(line, delim)
  unlist(words)
}
# =========================================================================== #
#' @title require_namespace
#' @description a wrapper for `requireNamespace`. Allowing us to more
#' easily create mock test for this functionality.
#' @param package_name the name of the package.
#' @returns A logical TRUE or FALSE boolean depending on whether the
#' package was added to the search path.
#' @keywords internal
#' @noRd
require_namespace <- function(package_name) {
  requireNamespace(package_name, quietly = TRUE)
}

# =========================================================================== #
