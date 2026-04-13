#' @title Get available processors
#' @name get_available_processors
#' @description
#' Get the number of available cores
#' @examples
#'
#' get_available_processors()
#'
#' @return Integer
#' @export
get_available_processors <- function() {
    parallelly::availableCores()
}
