#' @title read_qiime2_metadata
#' @description
#' Read a \href{https://qiime2.org}{qiime2} .tsv table containing metadata.
#'
#' @param metadata file name, a qiime2 .tsv file containing metadata about your
#' analysis.
#'
#' @examples
#'
#' metadata <- read_qiime2_metadata(rdataset_example(
#'   "sample_metadata.tsv"
#' ))
#'
#' @return A data.frame containing metadata
#' @export
read_qiime2_metadata <- function(metadata) {
  if (!file.exists(metadata)) {
    .abort_nonexistant_file(metadata)
  }


  data.frame()
}
