#' @title save_dataset
#' @description The save_dataset function will save the
#' \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#' to file.
#'
#' @param data a
#'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#'
#' @param file a string containing the file name.
#' @examples
#'
#' data <- read_mothur(
#'   fasta = strollur_example("final.fasta.gz"),
#'   count = strollur_example("final.count_table.gz"),
#'   taxonomy = strollur_example("final.taxonomy.gz"),
#'   design = strollur_example("mouse.time.design"),
#'   otu_list = strollur_example("final.opti_mcc.list.gz"),
#'   dataset_name = "miseq_sop"
#' )
#'
#' save_dataset(data, "miseq_sop.rds")
#'
#' @return A file containing the `strollur` object
#' @export
save_dataset <- function(data, file) {
  if (class(data)[1] != "strollur") {
    .abort_incorrect_type("strollur", data)
  }

  xint_serialize_dobject(data)
  saveRDS(data, file = file)
  file
}
