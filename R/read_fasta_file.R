#' @title read_fasta_file
#' @description
#' Read a FASTA formatted sequence file
#' @param filename FASTA file name (required)
#' @return A list containing the FASTA sequence data
#' @export
read_fasta_file = function(filename) {

    directory <- dirname(filename)
    filename <- basename(filename)

    fastafile <- file.path(file.path(directory), filename)
    df <- microseq::readFasta(fastafile)

    # extract name from comments
    names <- unlist(lapply(
        df$Header,
        (function(x) {
            val <- regexpr("\\s", x)
            if (val != -1) {
                name_comment <- substring(x, c(1, val + 1),
                                          c(val - 1, nchar(x)))
                return(name_comment[1])
            }
            return(x)
        })
    ))

    data <- list(names = names, sequences = df$Sequence)

    return(data)
}
