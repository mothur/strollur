#' @title read_mothur_shared
#' @description
#' Read a mothur formatted
#' \href{https://mothur.org/wiki/shared_file/}{shared file}
#' @param shared file name (required)
#' @examples
#'
#' # You can add your otu assignments to the your data set using the following:
#'
#' # read mothur shared file into data.frame
#' otu_data <- read_mothur_shared(strollur_example("final.opti_mcc.shared"))
#'
#' # create a new empty `strollur` object
#' data <- new_dataset()
#'
#' # assign abundance only 'otu' bins parsed by sample
#' assign(data = data, table = otu_data, type = "bins", bin_type = "otu")
#'
#' @return A data.frame containing the sequence otu assignments
#' @export
read_mothur_shared <- function(shared) {
  if (!file.exists(shared)) {
    .abort_nonexistant_file(shared)
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

  data.frame(
    bin_names = otu_assignments,
    abundances = abundance,
    samples = sample
  )
}
