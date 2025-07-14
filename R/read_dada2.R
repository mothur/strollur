#' @title read_dada2
#' @description
#' This function will read a dada2 sequence table and create a 'sequence_data'
#' object. The dada2 sequence table is a 2D matrix containing the abundance
#' counts by sample for each ASV. The sample names are stored as row names and
#' the sequence nucleotide strings are stored as column names.
#'
#' To generate the dada2 sequence table you can follow
#' \href{https://benjjneb.github.io/dada2/tutorial.html}{this dada2 tutorial}.
#' @param sequence_table A dada2 sequence table
#' @param dataset_name A string containing a name for your dataset.
#' @examples
#' \dontrun{
#' # Follow the dada2 tutorial above to generate seqtab from your fastq files.
#'
#' dim(seqtab)
#' [1]  20 293
#'
#' dataset <- read_dada2(seqtab, "dada2")
#'}
#' @return A 'sequence_data' object
#' @export
read_dada2 <- function(sequence_table, dataset_name = "") {

    # generate sequence names
    num_seqs <- ncol(sequence_table)
    seq_names <- c(1:num_seqs)
    seq_names <- paste("seq_", seq_names, sep="")

    # create new sequence_data object
    dataset <- sequence_data$new(dataset_name)
    dataset$add_sequences(seq_names,
                          colnames(sequence_table),
                          rep("dada2", num_seqs))

    sample_names <- rownames(sequence_table)
    names <- c()
    samples <- c()
    abundances <- c()

    # create a count table
    for (i in 1:ncol(sequence_table)) {
        abunds <- sequence_table[ ,i]
        non_zero_samples <- which(abunds != 0)

        names <- c(names, rep(seq_names[i], length(non_zero_samples)))
        samples <- c(samples, sample_names[non_zero_samples])
        abundances <- c(abundances, abunds[non_zero_samples])
    }

    dataset$assign_sequence_abundance(names, abundances, samples)
    dataset
}
