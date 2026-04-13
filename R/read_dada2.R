#' @title Create a
#'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#'   from dada2 outputs
#' @name read_dada2
#' @rdname read_dada2
#' @description
#' This function reads a dada2 sequence table and creates a `strollur`
#' object. The dada2 sequence table is a 2D matrix containing the abundance
#' counts by sample for each ASV. The sample names are stored as row names and
#' the sequence nucleotide strings are stored as column names.
#'
#' To generate the dada2 sequence table from your own files you can follow
#' \href{https://benjjneb.github.io/dada2/tutorial.html}{this dada2 tutorial}.
#'
#' @param sequence_table A dada2 sequence table
#' @param dataset_name A string containing a name for your dataset.
#' @examples
#'
#' seqtab <- readRDS(strollur_example("dada2.rds"))
#' dim(seqtab)
#'
#' data <- read_dada2(sequence_table = seqtab, dataset_name = "dada2 example")
#'
#' @return A `strollur` object
#' @export
read_dada2 <- function(sequence_table, dataset_name = "") {
  # generate sequence names
  num_seqs <- ncol(sequence_table)
  seq_names <- c(1:num_seqs)
  seq_names <- paste("seq_", seq_names, sep = "")

  # create new dataset object
  data <- new_dataset(dataset_name)
  xdev_add_sequences(data, data.frame(
    sequence_names = seq_names,
    sequences = colnames(sequence_table),
    comments = rep("dada2", num_seqs)
  ))

  sample_names <- rownames(sequence_table)
  names <- c()
  samples <- c()
  abundances <- c()

  # create a count table
  index <- 1
  for (col in colnames(sequence_table)) {
    abunds <- sequence_table[, index]
    non_zero_samples <- which(abunds != 0)

    names <- c(names, rep(seq_names[index], length(non_zero_samples)))
    samples <- c(samples, sample_names[non_zero_samples])
    abundances <- c(abundances, abunds[non_zero_samples])
    index <- index + 1
  }

  xdev_assign_sequence_abundance(data, data.frame(
    sequence_names = names,
    abundances = abundances,
    samples = samples
  ))
  data
}
