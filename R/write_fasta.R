#' @title write_fasta
#' @description
#'
#' Write a \href{https://www.ncbi.nlm.nih.gov/genbank/fastaformat/}{FASTA}
#' formatted sequence file
#'
#' @param data A 'dataset' object
#' @param filename a string containing the name of the output file. Default =
#' 'dataset_name'.fasta
#'
#' @examples
#'
#' miseq <- miseq_sop_example()
#' write_fasta(miseq)
#'
#' @return name of FASTA file
#' @export
write_fasta <- function(data, filename = NULL) {
  # check type
  if (class(data)[1] != "dataset") {
    .abort_incorrect_type("dataset", data)
  }

  if (is.null(filename)) {
    filename <- name(data, "dataset")
    if (filename == "") {
      .abort_no_name()
    }
    filename <- paste0(filename, ".fasta")
  }

  sequence_names <- name(data, "sequences")

  # data contains sequences
  if (length(sequence_names) != 0) {
    sequences <- get_sequences(data)

    # make sure they aren't blank
    if (!any(sequences == "")) {
      df <- data.frame(
        Header = sequence_names,
        Sequence = sequences
      )
      writeFasta(df, out.file = filename, width = 80)

      return(filename)
    }
  }

  return("no_sequence_data")
}
