#' @title read_mothur_count
#' @description
#' Read a mothur formatted count file
#' @param filename count file name (required)
#' @examples
#'
#' # mothur count file
#' # Representative_Sequence     total   sample2	sample3	sample4
#' # seq1	1150	250	400	500
#' # seq2	115	25	40	50
#' # seq3	50	25	25	0
#' # seq4	4	0	0	4
#'
#' # returns
#' # id   sample abundance
#' # <char>  <char>     <int>
#' #  1:   seq1 sample2       250
#' #  2:   seq1 sample3       400
#' #  3:   seq1 sample4       500
#' #  4:   seq2 sample2        25
#' #  5:   seq2 sample3        40
#' #  6:   seq2 sample4        50
#' #  7:   seq3 sample2        25
#' #  8:   seq3 sample3        25
#' #  9:   seq4 sample4         4
#'
#' # read a count file with samples
#' sample_table <- read_mothur_count(rdataset_example("final.count_table"))
#'
#' # You can add your sequence abundance data to your dataset as follows:
#'
#' dataset <- sequence_data$new()
#' dataset$add_sequences(unique(sample_table$id))
#' dataset$assign_sequence_abundance(
#'   sample_table$id,
#'   sample_table$abundance,
#'   sample_table$sample
#' )
#' dataset
#'
#' @return data.table
#' @export
read_mothur_count <- function(filename) {
  directory <- dirname(filename)
  filename <- basename(filename)

  count_file <- file.path(file.path(directory), filename)
  file_conn <- file(count_file)
  file_data <- readLines(file_conn)
  num_lines <- length(file_data)
  close(file_conn)

  table_names <- c()
  table_samples <- c()
  table_abunds <- c()

  # check for sample info in header
  if (num_lines > 1) {
    comment <- regexpr("#", file_data[2])

    # compressed
    if (comment[1] != -1) {
      # extract sample names
      # line 2 looks like: "#2,sample2	3,sample3	1,sample4"
      # remove first '#'
      pieces <- strsplit(file_data[2], "#")[[1]]
      file_data[2] <- pieces[nzchar(pieces)]
      words <- split_white_space(file_data[2])
      num_seqs <- length(file_data) - 3

      samples <- c()
      for (i in seq_along(words)) {
        # parse sample name
        file_index <- split_at_char(words[i], ",")

        # save sample names
        samples <- c(samples, file_index[2])
      }

      # read compressed data lines
      for (i in 4:length(file_data)) {
        seq_line <- split_white_space(file_data[i])

        # add copy of name for each abund
        table_names <- c(table_names, rep(
          seq_line[1],
          length(seq_line) - 2
        ))
        # skip name and total
        for (j in 3:length(seq_line)) {
          # looks like 2,3 -> meaning sample 2 has abundance 3
          data <- split_at_char(seq_line[j], ",")

          # add sample name to table_samples
          table_samples <- c(table_samples, samples[as.integer(data[1])])
          # add sample abundance to table_abunds
          table_abunds <- c(table_abunds, as.integer(data[2]))
        }
      }
    } else {
      # uncompressed format
      # Representative_Sequence  total  sample2	sample3	sample4
      words <- split_white_space(file_data[1])

      num_seqs <- length(file_data) - 1
      has_sample_data <- TRUE
      samples <- c()

      # no sample data in file
      if (length(words) == 2) {
        has_sample_data <- FALSE
        table_names <- rep("", num_seqs)
        table_abunds <- rep(0, num_seqs)
      } else {
        for (i in 3:length(words)) {
          samples <- c(samples, words[i])
        }
      }

      # read uncompressed data
      for (i in 2:length(file_data)) {
        seq_line <- split_white_space(file_data[i])

        if (length(seq_line) >= 2) {
          name <- seq_line[1]

          if (has_sample_data) {
            for (j in 3:length(seq_line)) {
              abund <- as.integer(seq_line[j])
              if (abund != 0) {
                table_names <- c(table_names, name)
                table_samples <- c(table_samples, samples[j - 2])
                table_abunds <- c(table_abunds, abund)
              }
            }
          } else {
            table_names[i - 1] <- name
            table_abunds[i - 1] <- as.integer(seq_line[2])
          }
        }
      }

      if (!has_sample_data) {
        return(data.table(
          id = table_names,
          abundance = table_abunds
        ))
      }
    }
  }

  data.table(
    id = table_names,
    sample = table_samples,
    abundance = table_abunds
  )
}
