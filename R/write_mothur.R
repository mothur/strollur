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
#'  written to files.
#' @param compress boolean, Default = TRUE.
#'
#' @examples
#'
#' miseq <- miseq_sop_example()
#' files <- write_mothur(miseq)
#'
#' @return a vector of file names
#' @export
write_mothur <- function(data, dir_path = NULL, compress = TRUE, tags = NULL) {
  # check type
  if (class(data)[1] != "dataset") {
    abort_incorrect_type("dataset", data)
  }

  dataset_name <- data$get_dataset_name()

  if (dataset_name == "") {
    dataset_name <- paste0("rdataset.", as.integer(Sys.time()),
      collapse = ""
    )
  }

  # if no dir given, set to current working directory
  if (is.null(dir_path)) {
    dir_path <- paste0(getwd(), .Platform$file.sep, dataset_name,
      collapse = ""
    )
  }

  if (!dir.exists(dir_path)) {
    # create it
    dir.create(dir_path)
  }

  # Interpret the result
  if (file.access(dir_path, mode = 2) != 0) {
    abort_not_writable(dir_path)
  }

  outputs <- c()

  # TODO
  ht <- !(is.null(tags))
  wrote_design <- FALSE
  wrote_count <- FALSE

  # write sequence_data, fasta, count, taxonomy
  if (!ht || ("sequence_data" %in% tags)) {
    filename <- paste0(dataset_name, ".fasta", collapse = "")
    output <- write_fasta(data, file.path(dir_path, filename))

    # will return "no_sequence_data" if no fasta data
    if (output != "no_sequence_data") {
      outputs <- c(outputs, output)
    }

    filename <- paste0(dataset_name, ".taxonomy", collapse = "")
    output <- write_taxonomy(data, file.path(dir_path, filename))

    # will return "no_sequence_taxonomy" if no classification data
    if (output != "no_sequence_taxonomy") {
      outputs <- c(outputs, output)
    }

    filename <- paste0(dataset_name, ".count_table", collapse = "")
    output <- write_mothur_count(data, file.path(dir_path, filename))

    # will return "no_sequence_data" if no sequence abundance data
    if (output != "no_sequence_data") {
      outputs <- c(outputs, output)
      wrote_count <- TRUE
    }

    if (data$get_num_treatments() != 0) {
      filename <- paste0(dataset_name, ".design", collapse = "")
      output <- write_mothur_design(data, file.path(dir_path, filename))
      outputs <- c(outputs, output)
      wrote_design <- TRUE
    }
  }

  # write bin_data
  if (!ht || ("bin_data" %in% tags)) {
    output <- write_mothur_list(data, file.path(dir_path, dataset_name))
    outputs <- c(outputs, output)

    # write count file, if needed
    if ((length(output) != 0) && !wrote_count) {
      filename <- paste0(dataset_name, ".count_table", collapse = "")
      output <- write_mothur_count(data, file.path(dir_path, filename))
      outputs <- c(outputs, output)
    }

    output <- write_mothur_shared(data, file.path(dir_path, dataset_name))
    outputs <- c(outputs, output)

    output <- write_mothur_cons_taxonomy(data, file.path(
      dir_path,
      dataset_name
    ))
    outputs <- c(outputs, output)

    if ((data$get_num_treatments() != 0) && !wrote_design) {
      filename <- paste0(dataset_name, ".design", collapse = "")
      output <- write_mothur_design(data, file.path(dir_path, filename))
      outputs <- c(outputs, output)
      wrote_design <- TRUE
    }
  }

  if (!ht || ("metadata" %in% tags)) {
    metadata <- data$get_metadata()

    if (nrow(metadata) != 0) {
      filename <- file.path(
        dir_path,
        paste0(dataset_name, ".metadata", collapse = "")
      )
      readr::write_tsv(metadata, filename)
      outputs <- c(outputs, filename)
    }
  }

  if (!ht || ("alignment_report" %in% tags)) {
    report <- data$get_alignment_report()

    if (nrow(report) != 0) {
      filename <- file.path(
        dir_path,
        paste0(dataset_name, ".alignment_report",
          collapse = ""
        )
      )
      readr::write_tsv(report, filename)
      outputs <- c(outputs, filename)
    }
  }

  if (!ht || ("contigs_assembly_report" %in% tags)) {
    report <- data$get_contigs_assembly_report()

    if (nrow(report) != 0) {
      filename <- file.path(
        dir_path,
        paste0(dataset_name, ".contigs_report",
          collapse = ""
        )
      )
      readr::write_tsv(report, filename)
      outputs <- c(outputs, filename)
    }
  }

  if (!ht || ("chimera_report" %in% tags)) {
    report <- data$get_chimera_report()

    if (nrow(report) != 0) {
      filename <- file.path(
        dir_path,
        paste0(dataset_name, ".chimera_report",
          collapse = ""
        )
      )
      readr::write_tsv(report, filename)
      outputs <- c(outputs, filename)
    }
  }

  if (!ht || ("references" %in% tags)) {
    references <- data$get_references()

    if (nrow(references) != 0) {
      filename <- file.path(
        dir_path,
        paste0(dataset_name, ".references",
          collapse = ""
        )
      )
      readr::write_tsv(references, filename)
      outputs <- c(outputs, filename)
    }
  }

  if (!ht || ("sample_tree" %in% tags)) {
    sample_tree <- data$get_sample_tree()

    if (!is.null(sample_tree)) {
      filename <- file.path(
        dir_path,
        paste0(dataset_name, ".sample.tree",
          collapse = ""
        )
      )
      ape::write.tree(sample_tree, file = filename)
      outputs <- c(outputs, filename)
    }
  }

  if (!ht || ("sequence_tree" %in% tags)) {
    sequence_tree <- data$get_sequence_tree()

    if (!is.null(sequence_tree)) {
      filename <- file.path(
        dir_path,
        paste0(dataset_name, ".sequence.tree", collapse = "")
      )
      ape::write.tree(sequence_tree, file = filename)
      outputs <- c(outputs, filename)
    }
  }

  if (compress) {
    # compress the files
    filename <- paste0(dir_path, ".zip", collapse = "")
    utils::zip(zipfile = filename, files = dir_path)
    unlink(dir_path, recursive = TRUE)
  }

  outputs
}
