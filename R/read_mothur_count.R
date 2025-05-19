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
#' # id   group abundance
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
#' # read a count file with groups
#' sample_table <- read_mothur_count(rdataset_example("test.count_table"))
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
  table_groups <- c()
  table_abunds <- c()

  # check for group info in header
  if (num_lines > 1) {
    comment <- regexpr("#", file_data[2])

    # compressed
    if (comment[1] != -1) {
      # extract group names
      # line 2 looks like: "#2,sample2	3,sample3	1,sample4"
      # remove first '#'
      pieces <- strsplit(file_data[2], "#")[[1]]
      file_data[2] <- pieces[nzchar(pieces)]
      words <- split_white_space(file_data[2])
      num_seqs <- length(file_data) - 3

      groups <- c()
      for (i in seq_along(words)) {
        # parse group name
        file_index <- split_at_char(words[i], ",")

        # save group names
        groups <- c(groups, file_index[2])
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
          # looks like 2,3 -> meaning group 2 has abundance 3
          data <- split_at_char(seq_line[j], ",")

          # add sample name to table_groups
          table_groups <- c(table_groups, groups[as.integer(data[1])])
          # add sample abundance to table_abunds
          table_abunds <- c(table_abunds, as.integer(data[2]))
        }
      }
    } else {
      # uncompressed format
      # Representative_Sequence  total  sample2	sample3	sample4
      words <- split_white_space(file_data[1])

      num_seqs <- length(file_data) - 1
      has_group_data <- TRUE
      groups <- c()

      # no group data in file
      if (length(words) == 2) {
        has_group_data <- FALSE
        table_names <- rep("", num_seqs)
        table_abunds <- rep(0, num_seqs)
      } else {
        for (i in 3:length(words)) {
          groups <- c(groups, words[i])
        }
      }

      # read uncompressed data
      for (i in 2:length(file_data)) {
        seq_line <- split_white_space(file_data[i])

        if (length(seq_line) >= 2) {
          name <- seq_line[1]

          if (has_group_data) {
            for (j in 3:length(seq_line)) {
              abund <- as.integer(seq_line[j])
              if (abund != 0) {
                table_names <- c(table_names, name)
                table_groups <- c(table_groups, groups[j - 2])
                table_abunds <- c(table_abunds, abund)
              }
            }
          } else {
            table_names[i - 1] <- name
            table_abunds[i - 1] <- as.integer(seq_line[2])
          }
        }
      }

      if (!has_group_data) {
        return(data.table(
          id = table_names,
          abundance = table_abunds
        ))
      }
    }
  }

  data.table(
    id = table_names,
    group = table_groups,
    abundance = table_abunds
  )
}
