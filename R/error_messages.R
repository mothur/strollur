#' @title abort_nonexistant_file
#' @description
#' Report file does not exist
#' @param filename String, containing name of the missing file
abort_nonexistant_file <- function(filename) {
  cli::cli_abort("[ERROR]: {.var {filename}} does not exist.")
}

#' @title abort_incorrect_type
#' @description
#' Report incorrect type provided
#' @param type_expected String, containing name of the type expected
#' @param type_received String, containing name of the type received
abort_incorrect_type <- function(type_expected, type_received) {
  cli::cli_abort("[ERROR]: Expected a {.var {type_expected}}
          but received a {.var {type_received}}.")
}

#' @title abort_provide_at_least_one
#' @description
#' Report missing parameter
#' @param parameters a list of parameter where one must be selected
abort_provide_at_least_one <- function(parameters) {
  message <- paste("[ERROR]: Expected at least one of the following parameters",
    "to be provided ", paste(parameters, collapse = ", "),
    ".",
    collapse = ""
  )
  cli::cli_abort(message)
}
