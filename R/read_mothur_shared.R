#' @title read_mothur_shared
#' @description
#' Read a mothur formatted
#' \href{https://mothur.org/wiki/shared_file/}{shared file}
#' @param shared file name (required)
#' @examples
#'
#' # You can add your otu assignments to the your data set using the following:
#'
#' otu_data <- read_mothur_shared(rdataset_example("final.opti_mcc.shared"))
#'
#' dataset <- sequence_data$new()
#' dataset$assign_bins(otu_data$bin_id, otu_data$abundance, otu_data$sample)
#'
#' @return A data.table containing the sequence otu assignments
#' @export
read_mothur_shared <- function(shared) {
  if (!file.exists(shared)) {
    abort_nonexistant_file(shared)
  }

  df <- readr::read_table(
    file = shared, col_names = TRUE,
    show_col_types = FALSE
  )

  sample_names <- df[[2]]
  # remove label, group and numOtus columns
  df <- df[, -c(1, 2, 3)]
  otu_names <- names(df)

  otu_assignments <- c()
  sample <- c()
  abundance <- c()

  i <- 1
  for (otu in df) {
    # only store non zero abundances
    non_zero_index <- which(otu != 0)
    sample <- c(sample, sample_names[non_zero_index])
    abundance <- c(abundance, otu[non_zero_index])
    otu_assignments <- c(
      otu_assignments,
      rep(otu_names[i], length(non_zero_index))
    )
    i <- i + 1
  }

  data.table(
    bin_id = otu_assignments,
    abundance = abundance,
    sample = sample
  )
}
