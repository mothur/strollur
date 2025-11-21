#' @title .abort_nonexistant_file
#' @description
#' Report directory / file does not exist
#' @param filename String, containing name of the missing file or directory
#' @keywords internal
.abort_nonexistant_file <- function(filename) {
  cli::cli_abort("[ERROR]: {.var {filename}} does not exist.")
}

#' @title .abort_no_name
#' @description
#' Report dataset has no name
#' @keywords internal
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
.abort_not_writable <- function(filename) {
  cli::cli_abort("[ERROR]: {.var {filename}} is not writable.")
}

#' @title .abort_incorrect_type
#' @description
#' Report incorrect type provided
#' @param type_expected String, containing name of the type expected
#' @param obj R object of wrong type
#' @keywords internal
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
.abort_incorrect_format <- function(type_expected, type_received) {
  cli::cli_abort("[ERROR]: Expected a {.var {type_expected}}
          but received {.var {type_received}}.")
}

#' @title .abort_provide_at_least_one
#' @description
#' Report missing parameter
#' @param parameters a list of parameter where one must be selected
#' @keywords internal
.abort_provide_at_least_one <- function(parameters) {
  message <- paste("[ERROR]: Expected at least one of the following parameters",
    "to be provided ", paste(parameters, collapse = ", "),
    ".",
    collapse = ""
  )
  cli::cli_abort(message)
}

#' @title abort_missing_column
#' @description
#' Data.frame is missing necessary columns
#' @param parameter name of missing required column
#' @export
abort_missing_column <- function(parameter) {
  message <- paste("Expected a data.frame column named ",
    parameter, " to be provided.",
    collapse = ""
  )
  cli::cli_abort(message)
}

#' @title alert_missing_data
#' @description
#' Import table is missing requested tag
#' @param tag name of missing requested tag
#' @keywords internal
.abort_missing_tag_alert <- function(tag) {
  message <- paste0(
    "[WARNING]: The import table does not include ",
    "'{.var {tag}', ignoring tag."
  )
  cli::cli_alert(message)
}

#' @title added_message
#' @description
#' Report dataset additions
#' @param num integer containing the number of items added.
#' @param tag string containing item added. Default = 'sequences'.
#' @export
added_message <- function(num = NULL, tag = "sequences") {
  if (!is.null(num)) {
    message <- paste0("Added ", as.character(num), " ", tag, ".",
      collapse = ""
    )
    cli::cli_alert_info(message)
  } else {
    message <- paste0("Added ", tag, ".", collapse = "")
    cli::cli_alert_info(message)
  }
}

#' @title assigned_message
#' @description
#' Report dataset assignments
#' @param num integer containing the number of items assigned
#' @param tag string containing item assigned Default = 'sequences'.
#' @export
assigned_message <- function(num, tag = "sequences") {
  message <- paste0("Assigned ", as.character(num), tag,
    collapse = ""
  )
  cli::cli_alert_info(message)
}
