#' @importFrom microseq readFastq
#' @title read_fastq
#'
#' @description
#' Read a \href{https://pmc.ncbi.nlm.nih.gov/articles/PMC2847217/}{FASTQ}
#' formatted file
#'
#' @param fastq string, containing the FASTQ file name. Accepts compressed `.gz`
#'   files as well.
#' @param format string, containing the FASTQ format. Options include:
#'   `illumina1.8+`, `solexa`, `sanger`. Default = `illumina1.8+`.
#' @param path string, optional path to FASTQ file.
#'
#' @examples
#'
#' fastq_data <- strollur::read_fastq(strollur_example("tiny.fastq.gz"))
#'
#' # fastq_data is a data.frame.
#' # To access the names of the sequences in the file, run the following:
#'
#' fastq_data$sequence_name
#'
#' # To access the sequences in the file, run the following:
#'
#' fastq_data$sequence
#'
#' # To access the quality scores in the file, run the following:
#'
#' fastq_data$quality_score
#'
#' @return A data.frame containing the FASTQ sequence data
#' @export
read_fastq <- function(fastq, path = NULL, format = "illumina1.8+") {
  if (!file.exists(fastq)) {
    .abort_nonexistant_file(fastq)
  }

  # if path is set, adjust fastq file name
  if (!is.null(path)) {
    fastq <- file.path(path, basename(fastq))
  }

  # use microseq to read fasta file
  df <- microseq::readFastq(fastq)

  # convert quality string to quality scores
  table <- NULL
  if (format == "solexa") {
    table <- .solexa_convert_table()$to_score
  }

  scores <- vector("list", length = length(df$Sequence))

  for (i in seq_along(df$Quality)) {
    scores[[i]] <- convert_qual(df$Quality[i], format, table)
  }

  data.frame(
    sequence_name = gsub(" .*$", "", df$Header),
    sequence = df$Sequence,
    quality_score = I(scores)
  )
}

#' @title .solexa_convert_table
#' @description Create table for solexa conversion
#'
#' @returns data.frame containing 2 columns to_score and to_string for
#'   converting solexa quality string to sanger quality strings
#'
#' @keywords internal
#' @noRd
.solexa_convert_table <- function() {
  i <- -64:64
  raw_values <- 33 + 10 * log10(1 + 10^(i / 10.0)) + 0.499
  int_to_char <- as.integer(raw_values)
  char_to_int <- intToUtf8(int_to_char, multiple = TRUE)
  data.frame(to_score = char_to_int, to_string = int_to_char)
}

#' @title convert_qual
#' @description Convert quality string to quality scores
#'
#' @param quality_string, string containing quality data string
#' @param format string, containing the FASTQ format. Options include:
#'   `illumina1.8+`, `solexa`, `sanger`. Default = `illumina1.8+`.
#' @param table quality conversion table for solexa format. Default = `NULL`.
#' @examples
#'
#' quality_string <- "3AA?ABBDBFFBEGGEGGGGAFFGGGGGHHHCGGGGGGHFGHGGCFDEFGGGH"
#' convert_qual(quality_string)
#'
#' @returns vector of integers
#' @export
convert_qual <- function(quality_string,
                         format = "illumina1.8+",
                         table = NULL) {
  quality_scores <- utf8ToInt(quality_string)
  has_negative_scores <- FALSE

  if (is.null(table) && (format == "solexa")) {
    table <- .solexa_convert_table()$to_score
  }

  for (i in seq_along(quality_scores)) {
    temp <- quality_scores[i]
    if (format == "illumina") {
      temp <- temp - 64 # char '@'
    } else if (format == "illumina1.8+") {
      temp <- temp - 33 # char '!' //33
    } else if (format == "solexa") {
      temp <- utf8ToInt(table[temp]) # convert to sanger
      temp <- temp - 33 # char '!' //33
    } else {
      temp <- temp - 33 # char '!' //33
    }

    if (temp < 0) {
      has_negative_scores <- TRUE
      temp <- 0
    }
    quality_scores[i] <- temp
  }

  if (has_negative_scores) {
    message <- paste0(
      "Found negative quality scores, do you have the ",
      "right format selected? Please refer to FASTQ ",
      "documentation."
    )
    cli::cli_abort(message)
  }

  quality_scores
}
