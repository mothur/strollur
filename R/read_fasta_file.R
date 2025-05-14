#' @title read_fasta_file
#' @description
#' Read a FASTA formatted sequence file
#' @param fasta FASTA file name (required)
#' @examples
#'
#' fasta_data <- read_fasta_file(rdataset_example("test.fasta"))
#'
#' # fasta_data is a named list.
#' # To access the names of the sequences in the file, run the following:
#'
#' fasta_data$names
#'
#' # To access the sequences in the file, run the following:
#'
#' fasta_data$sequences
#'
#' @return A list containing the FASTA sequence data
#' @export
read_fasta_file <- function(fasta) {
  if (!file.exists(fasta)) {
      abort_nonexistant_file(fasta)
  }

  # use microseq to read fasta file
  df <- microseq::readFasta(fasta)

  # extract name from comments
  num_seqs <- nrow(df)

  has_comments <- FALSE
  names <- rep("", num_seqs)
  comments <- rep("", num_seqs)

  for (i in 1:num_seqs) {
    names[i] <- df$Header[i]

    val <- regexpr("\\s", names[i])
    if (val != -1) {
      name_comment <- substring(
        names[i], c(1, val + 1),
        c(val - 1, nchar(names[i]))
      )
      names[i] <- name_comment[1]
      comments[i] <- name_comment[2]
      has_comments <- TRUE
    }
  }

  if (has_comments) {
    return(list(names = names, sequences = df$Sequence, comments = comments))
  }

  list(names = names, sequences = df$Sequence)
}
