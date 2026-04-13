# Define an S3 generic - this allows for the additional parameters to names
names <- function(x, type = "sequences",
                  bin_type = "otu",
                  sample = NULL,
                  distinct = FALSE) {
  UseMethod("names", x)
}

#' @title Get the names of various data
#'   in a \href{https://mothur.org/strollur/reference/strollur.html}{strollur}
#'   object
#' @name names
#' @rdname names
#' @description Get the names of names sequences, bins, samples, treatments, and
#'   reports data in a
#'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#'
#' @param data, a
#'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#'
#' @param type, string containing the type of data you would like. Options
#' include: "dataset", "sequences", "bins", "samples", "treatments", "reports".
#' Default = "sequences".
#'
#' @param bin_type, string containing the bin type you would like the names
#' for. Default = "otu".
#'
#' @param samples, vector of strings. samples is only used when 'type' =
#' "sequences" or 'type' = "bins" . samples should contain the names of the
#' samples you want names for. Default = NULL.
#'
#' @param distinct, Boolean. distinct is used when 'type' =
#' "sequences" or 'type' = "bins" and the samples parameter is used. The
#' distinct parameter allows you to get the names that present given
#' set of samples. When distinct is TRUE, the names function will return the
#' names that ONLY contain data from the given samples. When distinct is FALSE
#' the data returned contains data from a given samples, but may ALSO contain
#' data from other samples. Default = FALSE.
#'
#' @examples
#'
#' miseq <- miseq_sop_example()
#'
#' # To get the name of the dataset
#' names(data = miseq, type = "dataset")
#'
#' # To get the names of the sequences
#' names(data = miseq, type = "sequences")
#'
#' # To get the names of the sequences present sample 'F3D0'
#' names(data = miseq, type = "sequences", samples = c("F3D0"))
#'
#' #' # To get the names of the sequences unique to sample 'F3D0'
#' names(data = miseq, type = "sequences", samples = c("F3D0"), distinct = TRUE)
#'
#' # To get the names of the samples
#' names(data = miseq, type = "samples")
#'
#' # To get the names of the treatments
#' names(data = miseq, type = "treatments")
#'
#' # To get the names of the bins
#' names(data = miseq, type = "bins")
#'
#' # To get the names of the bins that are unique to 'F3D0'
#' names(data = miseq, type = "bins", samples = c("F3D0"), distinct = TRUE)
#'
#' # To get the names of the bins that include sequences from 'F3D0'
#' names(data = miseq, type = "bins", samples = c("F3D0"), distinct = FALSE)
#'
#' # To get the names of the reports
#' names(data = miseq, type = "reports")
#'
#' @return vector of strings, containing the names requested
#' @export
names <- function(data,
                  type = "sequences",
                  bin_type = "otu",
                  samples = NULL,
                  distinct = FALSE) {
  if (inherits(data, "strollur")) {
    # allow for type and bin_type to be entered without ""
    type <- as.character(substitute(type))
    bin_type <- as.character(substitute(bin_type))

    xdev_names(data, type, bin_type, samples, distinct)
  } else {
    base::names(data)
  }
}
