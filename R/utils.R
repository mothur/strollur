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
.split_at_char <- function(line, delim = ",") {
  words <- strsplit(line, delim)
  unlist(words)
}
# =========================================================================== #
