#' @title .abort_nonexistant_file
#' @description
#' Report directory / file does not exist
#' @param filename String, containing name of the missing file or directory
#' @keywords internal
#' @noRd
.abort_nonexistant_file <- function(filename) {
  cli::cli_abort("[ERROR]: {.var {filename}} does not exist.")
}

#' @title .abort_no_name
#' @description
#' Report dataset has no name
#' @keywords internal
#' @noRd
.abort_no_name <- function() {
  message <- paste0(
    "[ERROR]: You did not set an output file name, and your ",
    "dataset not not have a name associated with it. Please ",
    "provide an output file name."
  )
  cli::cli_abort(message)
}

#' @title .abort_not_writable
#' @description
#' Report directory / file is not writable
#' @param filename String, containing name of the file without write
#' permissions
#' @keywords internal
#' @noRd
.abort_not_writable <- function(filename) {
  cli::cli_abort("[ERROR]: {.var {filename}} is not writable.")
}

#' @title .abort_incorrect_type
#' @description
#' Report incorrect type provided
#' @param type_expected String, containing name of the type expected
#' @param obj R object of wrong type
#' @keywords internal
#' @noRd
.abort_incorrect_type <- function(type_expected, obj) {
  cli::cli_abort("[ERROR]: Expected a {.var {type_expected}}
          but received {.cls { class(obj)} }.")
}

#' @title .abort_incorrect_format
#' @description
#' Report incorrect artifact format provided
#' @param type_expected String, containing name of the type expected
#' @param type_received String, containing name of the type received
#' @keywords internal
#' @noRd
.abort_incorrect_format <- function(type_expected, type_received) {
  cli::cli_abort("[ERROR]: Expected a {.var {type_expected}}
          but received {.var {type_received}}.")
}
