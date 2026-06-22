#' @title read_mothur_oligos
#' @description
#' Read a mothur formatted
#' \href{https://mothur.org/wiki/oligos_file/}{oligos file}
#' @param oligos file name. a mothur formatted
#' \href{https://mothur.org/wiki/oligos_file/}{oligos file}
#' @examples
#'
#' oligos <- read_mothur_oligos(strollur_example("paired_read.oligos"))
#'
#' # Create a new dataset and add your oligos data
#'
#' data <- new_dataset() |>
#'   add(
#'     table = oligos,
#'     type = "report",
#'     report_type = "paired_oligos"
#'   )
#'
#' @return A data.frame containing the oligos data.
#' @export
read_mothur_oligos <- function(oligos) {
  if (!file.exists(oligos)) {
    .abort_nonexistant_file(oligos)
  }

  # read file
  oligos_lines <- readLines(oligos)

  # linkers (2 columns) : 'linker' linkerString
  # spacers (2 columns) : 'spacer' spacerString
  # forward (2 columns) : 'forward' primerString optionalName
  # reverse (2 columns) : 'reverse' primerString
  # barcode (3 columns) : 'barcode' barcodeString barcodeName
  # paired primers (4 column) : 'primer' forPrimer revPrimer primerName
  #                     : primerName is optional, NONE keyword is allowed
  # paired barcodes (4 column) : 'barcode' forBarcode revBarcode barcodeName
  #                     : barcodeName is required, NONE keyword is allowed

  # error if mixing paired and unpaired oligos types

  num_oligo <- length(oligos_lines)

  linker <- spacer <- forward <- reverse <- barcode <- barcode_name <- c()
  forward_name <- primer_forward <- primer_reverse <- c()
  primer_name <- barcode_forward <- barcode_reverse <- c()

  has_paired <- FALSE
  has_single <- FALSE
  for (oligo in oligos_lines) {
    oligo <- .split_white_space(oligo)
    oligo_length <- length(oligo)

    if (oligo[1] == "barcode") {
      if (oligo_length == 4) {
        has_paired <- TRUE
        barcode_forward <- c(barcode_forward, oligo[2])
        barcode_reverse <- c(barcode_reverse, oligo[3])
        barcode_name <- c(barcode_name, oligo[4])
      } else if (oligo_length == 3) {
        has_single <- TRUE
        barcode <- c(barcode, oligo[2])
        barcode_name <- c(barcode_name, oligo[3])
      } else {
        cli::cli_abort("barcode names are required")
      }
    } else if (oligo[1] == "primer") {
      if (oligo_length == 3) {
        has_paired <- TRUE
        primer_forward <- c(primer_forward, oligo[2])
        primer_reverse <- c(primer_reverse, oligo[3])
      } else if (oligo_length == 4) {
        has_paired <- TRUE
        primer_forward <- c(primer_forward, oligo[2])
        primer_reverse <- c(primer_reverse, oligo[3])
        primer_name <- c(primer_name, oligo[4])
      } else {
        cli::cli_abort(paste(
          "primers must be paired. You can use",
          "'NONE' as a placeholder for paired reads",
          "where only the forward or reverse primer",
          "is present. Use the `forward` or",
          "`reverse` type for unpaired primers."
        ))
      }
    } else if (oligo[1] == "forward") {
      has_single <- TRUE
      forward <- c(forward, oligo[2])
      forward_name <- c(
        forward_name,
        ifelse(length(oligo) == 3, oligo[3], "")
      )
    } else if (oligo[1] == "reverse") {
      has_single <- TRUE
      reverse <- c(reverse, oligo[2])
    } else if (oligo[1] == "linker") {
      has_single <- TRUE
      linker <- c(linker, oligo[2])
    } else if (oligo[1] == "spacer") {
      has_single <- TRUE
      spacer <- c(spacer, oligo[2])
    } else {
      cli::cli_abort(paste0(
        "{.var {oligo[1]}} is not a valid oligo ",
        "type. Options include: linker, spacer, ",
        "forward, reverse, barcode and primer."
      ))
    }
  }

  if (has_single && has_paired) {
    cli::cli_abort(paste0("cannot mix paired and unpaired oligo data"))
  } else if (has_single) {
    oligo_list <- list(
      barcode = barcode,
      barcode_name = barcode_name,
      forward = forward,
      forward_name = forward_name,
      reverse = reverse,
      linker = linker,
      spacer = spacer
    )
    df <- .convert_to_dataframe(oligo_list)
  } else if (has_paired) {
    oligo_list <- list(
      barcode_forward = barcode_forward,
      barcode_reverse = barcode_reverse,
      barcode_name = barcode_name,
      primer_forward = primer_forward,
      primer_reverse = primer_reverse,
      primer_name = primer_name,
      linker = linker,
      spacer = spacer
    )
    df <- .convert_to_dataframe(oligo_list)
  } else {
    df <- data.frame()
  }
  df
}

#' @keywords internal
#' @noRd
.convert_to_dataframe <- function(oligos_list) {
  # remove empty categories
  oligos_list <- oligos_list[sapply(oligos_list, length) > 0]
  # find oligos with longest list
  max_len <- max(sapply(oligos_list, length))
  # add blank lines to even out categories and create data.frame
  df <- as.data.frame(lapply(
    oligos_list,
    function(x) c(x, rep("", max_len - length(x)))
  ))
}
