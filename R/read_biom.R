#' @title Create a
#'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#'   from a biom formatted file.
#' @name read_biom
#' @rdname read_biom
#' @description
#' The read_biom function reads a \href{https://biom-format.org}{biom formatted}
#' file and creates a `strollur::strollur` object.
#'
#' @param biom string containing the name of the biom file.
#'
#' @examples
#'
#' if (requireNamespace("h5lite", quietly = TRUE)) {
#'   data <- strollur::read_biom(biom = strollur_example("miseq_sop.otu.biom"))
#'   data
#' } else {
#'   message(paste(
#'     "To use this functionality you have to install the",
#'     "h5lite package."
#'   ))
#' }
#'
#' @return A
#'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#' @export
read_biom <- function(biom) {
  if (!require_namespace("h5lite")) {
    stop(paste(
      "The h5lite R package is required to read HDF5 formatted",
      "BIOM files. To install h5lite, run: pak::pak('h5lite')."
    ))
  }

  if (!file.exists(biom)) {
    .abort_nonexistant_file(biom)
  }

  # read biom file
  hdata <- rbiom::read_biom(biom)

  data <- strollur::new_dataset(hdata$id)

  counts_matrix <- hdata$counts
  sample_indices <- rep(
    seq_len(ncol(counts_matrix)),
    diff(counts_matrix@p)
  )

  bin_names <- unique(hdata$otus[counts_matrix@i + 1])

  if (!is.null(hdata$sequences)) {
    add(data,
      table = data.frame(
        bin_name = bin_names,
        sequence = as.vector(hdata$sequences)
      ),
      type = "sequence",
      table_names = list(sequence_name = "bin_name")
    )
  } else {
    add(data,
      table = data.frame(bin_name = bin_names),
      type = "sequence",
      table_names = list(sequence_name = "bin_name")
    )
  }

  table <- data.frame(
    bin_name  = hdata$otus[counts_matrix@i + 1],
    abundance = counts_matrix@x,
    sample    = hdata$samples[sample_indices]
  )

  # assign sample and bin abundance - shared data
  xdev_assign_sequence_abundance(
    data = data,
    table = table,
    sequence_name = "bin_name"
  )

  # assign to bins
  xdev_assign_bins(data, table = table, bin_type = "asv")

  if (!is.null(hdata$taxonomy)) {
    # convert from factors to characters
    tax_character <- data.frame(
      lapply(hdata$taxonomy, as.character),
      stringsAsFactors = FALSE
    )

    # skip otu names in column 1
    taxes <- paste0(
      apply(tax_character[, 2:7], 1, paste, collapse = ";"),
      ";"
    )

    xdev_assign_bin_taxonomy(data,
      table = data.frame(
        bin_name = bin_names,
        taxonomy = taxes
      ),
      bin_type = "asv"
    )
  }

  if (!is.null(hdata$metadata)) {
    xdev_add_report(data, table = hdata$metadata)
  }

  data
}
