# Define an S3 generic - this allows for the additional parameters to count
count <- function(x, type = "sequences",
                  bin_type = "otu",
                  sample = NULL,
                  distinct = FALSE) {
  UseMethod("count", x)
}

#' @title count
#' @description
#' Find the number of sequences, samples, treatments or bins of a given type in
#' a \link{strollur} object
#'
#' @param data, a \link{strollur} object
#'
#' @param type, string containing the type of data you want the number of.
#' Options include: "sequences", "samples", "treatments", "bins", and
#'  "references". Default = "sequences".
#'
#' @param bin_type, string containing the bin type you would like the number of
#' bins for. Default = "otu".
#'
#' @param samples, vector of strings. samples is only used when 'type' =
#' "sequences" or 'type' = "bins" . samples should contain the names of the
#' samples you want the count for. Default = NULL.
#'
#' @param distinct, Boolean. distinct is used when 'type' =
#' "sequences" or 'type' = "bins". When 'type' = "sequences" and distinct is
#' TRUE the number of unique sequences is returned. When 'type' = "sequences"
#' and distinct is FALSE the total number of sequences is returned. This can
#' also be combined with samples to find the number of unique sequences found
#' ONLY in a given set of samples, or to find the number of unique sequences
#' in given set of samples that may also be present in other samples.
#' When 'type' = "bins", you can set distinct = TRUE to return the number of
#' bins that ONLY contain sequences from the given samples. When distinct is
#' FALSE the count returned contains bins with sequences from a given samples,
#' but those bins may also contain other samples. Default = FALSE.
#'
#' @examples
#'
#' miseq <- miseq_sop_example()
#'
#' # To get the total number of sequences
#' count(data = miseq, type = "sequences")
#'
#' # To get number of unique sequences
#' count(data = miseq, type = "sequences", distinct = TRUE)
#'
#' # To get number of unique sequences from samples 'F3D0' and 'F3D1'
#' # Note these sequences will be present in both samples but may be
#' # be present in other samples as well
#' count(data = miseq, type = "sequences", samples = c("F3D0", "F3D1"))
#'
#' # To get number of unique sequences exclusive to samples 'F3D0' and 'F3D1'
#' # Note sequences are present in both samples and NOT present in any other
#' # samples.
#' count(
#'   data = miseq, type = "sequences", samples = c("F3D0", "F3D1"),
#'   distinct = TRUE
#' )
#'
#' # To get the number of samples in the dataset
#' count(data = miseq, type = "samples")
#'
#' # To get the number of treatments in the dataset
#' count(data = miseq, type = "treatments")
#'
#' # To get the number of "otu" bins in the dataset
#' count(data = miseq, type = "bins", bin_type = "otu")
#'
#' # To get the number of "asv" bins in the dataset
#' count(data = miseq, type = "bins", bin_type = "asv")
#'
#' # To get the number of "phylotype" bins in the dataset
#' count(data = miseq, type = "bins", bin_type = "phylotype")
#'
#' # To get number of "otu" bins from samples 'F3D0' and 'F3D1'
#' # Note these bins will have sequences from both samples but there may be
#' # other samples present as well
#' count(
#'   data = miseq,
#'   type = "bins", bin_type = "otu", samples = c("F3D0", "F3D1")
#' )
#'
#' # To get number of "otu" bins unique to samples 'F3D0' and 'F3D1'
#' # Note these bins will have sequences from both samples and NO other samples
#' # will be present in the bins.
#' count(
#'   data = miseq, type = "bins", bin_type = "otu",
#'   samples = c("F3D0", "F3D1"), distinct = TRUE
#' )
#'
#' @return double
#' @export
count <- function(data,
                  type = "sequences",
                  bin_type = "otu",
                  samples = NULL,
                  distinct = FALSE) {
  if ("strollur" %in% class(data)) {
    xdev_count(data, type, bin_type, samples, distinct)
  } else {
    dplyr::count(data)
  }
}
