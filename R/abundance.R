#' @title Get the abundance data for sequences, bins, samples, and treatments in
#'   a \link{strollur::strollur} object
#' @name abundance
#' @rdname abundance
#' @description
#' Get the abundance data for sequences, bins, samples, and treatments in a
#' \link{strollur::strollur} object
#'
#' @param data, a \link{strollur::strollur} object
#'
#' @param type, string containing the type of data you want the number of.
#' Options include: "sequences", "bins", "samples" and "treatments".
#' Default = "sequences".
#'
#' @param bin_type, string containing the bin type you would like the abundance
#' data for. Default = "otu".
#'
#' @param by_sample, Boolean. When by_sample is TRUE, the abundance data will
#' be parsed by sample. Default = FALSE.
#'
#' @examples
#'
#' miseq <- miseq_sop_example()
#'
#' # To the total abundance for each sequence
#' abundance(data = miseq, type = "sequences")
#'
#' # To the total abundance for each sequence parsed by sample
#' abundance(data = miseq, type = "sequences", by_sample = TRUE)
#'
#' # To the total abundance for each "otu" bin
#' abundance(data = miseq, type = "bins", bin_type = "otu")
#'
#' # To the total abundance for each "otu" bin parsed by sample
#' abundance(data = miseq, type = "bins", bin_type = "otu", by_sample = TRUE)
#'
#' # To the total abundance for each "asv" bin
#' abundance(data = miseq, type = "bins", bin_type = "asv")
#'
#' # To the total abundance for each "asv" bin parsed by sample
#' abundance(data = miseq, type = "bins", bin_type = "asv", by_sample = TRUE)
#'
#' # To the total abundance for each sample
#' abundance(data = miseq, type = "samples")
#'
#' # To the total abundance for each treatment
#' abundance(data = miseq, type = "treatments")
#'
#' @return data.frame
#' @export
abundance <- function(data,
                      type = "sequences",
                      bin_type = "otu",
                      by_sample = FALSE) {
  xdev_abundance(data, type, bin_type, by_sample)
}
