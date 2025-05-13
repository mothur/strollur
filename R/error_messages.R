#' @title abort_name_missing
#' @description
#' Report name missing and abort
#' @param name String, containing name of missing sequence
abort_name_missing = function(name) {
    cli::cli_abort("[ERROR]: {.var {name}} is not in your dataset,
                            object mismatch.")
}

#' @title abort_length_mismatch
#' @description
#' Report length mismatch and abort
#' @param name1 String, containing name of first variable
#' @param name2 String, containing name of second variable
#' @param length1 Numeric, containing length of first variable
#' @param length2 Numeric, containing length of second variable
abort_length_mismatch = function(name1, name2, length1, length2) {
    cli::cli_abort("[ERROR]: The length of {.var {name1}}
          must equal the length of {.var {name2}}. {.var {name1}} has length
          {.var {length1}} but {.var {name2}} has length {.var {length2}}")
}

#' @title alert_missing_name
#' @description
#' Report missing name
#' @param name String, containing name of missing sequence
#' @param return_value String, containing return_value
alert_missing_name = function(name, return_value) {
    cli::cli_alert("{.var {name}} is not in your dataset,
                       returning {.var {return_value}}.")
}

#' @title alert_eliminated_name
#' @description
#' Report eliminated name
#' @param name String, containing name of eliminated sequence
#' @param return_value String, containing return_value
alert_eliminated_name = function(name, return_value) {
    cli::cli_alert("{.var {name}} has been eliminated from your
          dataset, returning {.var {return_value}}.")
}
