#' @title write_fasta
#' @description
#'
#' Write a \href{https://www.ncbi.nlm.nih.gov/genbank/fastaformat/}{FASTA}
#' formatted sequence file
#'
#' @param data A `strollur` object
#' @param filename a string containing the name of the output file. Default =
#' 'dataset_name'.fasta
#' @param degap a logical. Default = FALSE. When degap = `TRUE`, all gap
#' characters will be removed from the sequences.
#'
#' @examples
#'
#' miseq <- miseq_sop_example()
#' write_fasta(miseq, tempfile())
#'
#' @return name of FASTA file
#' @export
write_fasta <- function(data, filename = NULL, degap = FALSE) {
  if (!inherits(data, "strollur")) {
    stop("data must be a strollur object.")
  }

  if (is.null(filename)) {
    filename <- names(data, "dataset")
    if (filename == "") {
      .abort_no_name()
    }
    filename <- paste0(filename, ".fasta")
  }

  sequence_names <- names(data, "sequences")

  # data contains sequences
  if (length(sequence_names) != 0) {
    sequences <- xdev_get_sequences(data, degap = degap)

    # make sure they aren't blank
    if (!any(sequences == "")) {
      df <- data.frame(
        Header = sequence_names,
        Sequence = sequences
      )
      microseq::writeFasta(df, out.file = filename, width = 80)

      return(filename)
    }
  }

  "no_sequence_data"
}
