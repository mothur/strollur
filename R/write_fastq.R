#' @importFrom microseq writeFastq
#'
#' @title Write a
#'   \href{https://pmc.ncbi.nlm.nih.gov/articles/PMC2847217/}{FASTQ} formatted
#'   file
#' @description
#'
#' Write a \href{https://pmc.ncbi.nlm.nih.gov/articles/PMC2847217/}{FASTQ}
#' formatted file
#'
#' @param data a
#'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#' @param filename a string containing the name of the output file. Default =
#' 'dataset_name'.fastq
#' @examples
#'
#' table <- strollur::read_fastq(strollur_example("tiny.fastq.gz"))
#'
#' # Create `strollur::strollur` object
#' data <- strollur::new_dataset("example")
#'
#' # Add FASTQ data `strollur::strollur` object
#' strollur::add(data, table, type = "fastq")
#'
#' # Write `strollur::strollur` objects FASTQ data to file
#' strollur::write_fastq(data, tempfile())
#'
#' @return name of FASTQ file
#' @export
write_fastq <- function(data, filename = NULL) {
  if (!inherits(data, "strollur")) {
    stop("data must be a strollur object.")
  }

  if (is.null(filename)) {
    filename <- names(data, "dataset")
    if (filename == "") {
      .abort_no_name()
    }
    filename <- paste0(filename, ".fastq")
  }

  fastq <- xdev_report(data, type = "fastq")

  # data contains sequences
  if (nrow(fastq) != 0) {
    # look for fastq format
    fastq_format <- attr(fastq, "fastq_format")
    if (is.null(fastq_format)) {
      fastq_format <- "illumina1.8+"
    }

    table <- NULL
    if (fastq_format == "solexa") {
      i <- -64:64
      raw_values <- 33 + 10 * log10(1 + 10^(i / 10.0)) + 0.499
      table <- as.integer(raw_values)
    }

    df <- data.frame(
      Header = fastq$sequence_name,
      Sequence = fastq$sequence,
      Quality = sapply(fastq$quality_score, convert_qual_string)
    )
    microseq::writeFastq(df, out.file = filename)

    return(filename)
  }

  "no_fastq_data"
}

#' @title convert_qual_string
#' @description Convert quality scores to quality string
#'
#' @param quality_score vector containing the quality scores
#' @param format string, containing the FASTQ format. Options include:
#'   `illumina1.8+`, `solexa`, `sanger`. Default = `illumina1.8+`.
#' @param table quality conversion table for solexa format. Default = `NULL`.
#' @examples
#'
#' quality_score <- c(
#'   2, 29, 29, 32, 32, 33, 32, 33, 33,
#'   37, 37, 37, 38, 38, 38
#' )
#'
#' convert_qual_string(quality_score)
#'
#' @returns vector of characters
#' @keywords internal
#' @noRd
convert_qual_string <- function(quality_score,
                                format = "illumina1.8+",
                                table = NULL) {
  if (is.null(table) && (format == "solexa")) {
    i <- -64:64
    raw_values <- 33 + 10 * log10(1 + 10^(i / 10.0)) + 0.499
    table <- as.integer(raw_values)
  }


  control_char <- if (format == "illumina") utf8ToInt("@") else utf8ToInt("!")
  temp <- quality_score + control_char

  if (format == "solexa") {
    temp <- table[temp]
  }

  qual_string <- intToUtf8(temp)
  qual_string
}
