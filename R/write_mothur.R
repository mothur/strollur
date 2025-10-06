#' @title write_mothur
#' @description
#' The write_mothur function will write various
#' \href{https://mothur.org/wiki/tags/#file_types}{file types} for use with
#' mothur.
#'
#' @param data A 'dataset' object
#' @param dir_path a string containing the name of directory where the files
#' should be written. Default = current working directory.
#' @param tags a vector of strings containing the items you wish to write
#' Options are 'sequence_data', 'bin_data', 'metadata',
#' 'references', 'sequence_tree', 'sample_tree', 'alignment_report',
#' 'contigs_assembly_report' and 'chimera_report'. By default, everything is
#'  written to file.
#' @examples
#'
#' miseq <- miseq_sop_example()
#' files <- write_mothur(miseq)
#'
#' @return a vector of file names
#' @export
write_mothur <- function(data, dir_path = NULL, tags = NULL) {
  # check type
  if (class(data)[1] != "dataset") {
    abort_incorrect_type("dataset", data)
  }

  # if no dir given, set to current working directory
  if (is.null(dir_path)) {
    dir_path <- getwd()
  } else {
    if (!dir.exists(dir_path)) {
      # create it
      dir.create(dir_path)
    }

    # Interpret the result
    if (file.access(dir_path, mode = 2) != 0) {
      abort_not_writable(dir_path)
    }
  }

  # TODO
  if (is.null(tags)) {

  } else {

  }
}
