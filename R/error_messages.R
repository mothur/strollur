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
#' @param obj R object of wrong type
abort_incorrect_type <- function(type_expected, obj) {
  cli::cli_abort("[ERROR]: Expected a {.var {type_expected}}
          but received {.cls { class(obj)} }.")
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

#' @title abort_missing_column
#' @description
#' Data.frame is missing necessary columns
#' @param parameter name of missing required column
abort_missing_column <- function(parameter) {
  message <- paste("[ERROR]: Expected a data.frame column named ",
    parameter, " to be provided.",
    collapse = ""
  )
  cli::cli_abort(message)
}

#' @title added_message
#' @description
#' Report dataset additions
#' @param num integer containing the number of items added.
#' @param tag string containing item added. Default = 'sequences'.
added_message <- function(num, tag = "sequences") {
    message <- paste0("Added ", as.character(num), " ", tag, ".",
                     collapse = ""
    )
    cli::cli_alert_info(message)
}

#' @title assigned_message
#' @description
#' Report dataset assignments
#' @param num integer containing the number of items assigned
#' @param tag string containing item assigned Default = 'sequences'.
assigned_message <- function(num, tag = "sequences") {
    message <- paste0("Assigned ", as.character(num), tag,
                     collapse = ""
    )
    cli::cli_alert_info(message)
}
