#' @title read_qiime2_taxonomy
#' @description
#' Read a \href{https://qiime2.org}{qiime2} qza containing taxonomy data
#'
#' # nolint start
#' To generate the various input files you can follow \href{https://amplicon-docs.qiime2.org/en/latest/tutorials/moving-pictures.html}{qiime moving-pictures}.
#' # nolint end
#'
#' @param qza file name, a qiime2 .qza file containing taxonomy data.
#' @param dir_path a string containing the name of directory where the artifacts
#' files should be unpacked. Default = current working directory.
#' @param remove_unpacked_artifacts boolean, When TRUE, the artifact's
#' temporary directories will be removed after processing. Default = TRUE.
#'
#' @examples
#'
#' artifact <- read_qiime2_taxonomy(rdataset_example(
#'   "taxonomy.qza"
#' ))
#'
#' # access the taxonomy table
#'
#' artifact$data
#'
#' @return A list containing artifact
#' @export
read_qiime2_taxonomy <- function(qza, dir_path = NULL,
                                 remove_unpacked_artifacts = TRUE) {
  # if no dir_path given, set to current working directory
  if (is.null(dir_path)) {
    dir_path <- getwd()
  }

  if (!dir.exists(dir_path)) {
    # create it
    dir.create(dir_path)
  }

  artifact_name <- sub("\\.[^.]+$", "", basename(qza))
  tmp_path <- file.path(dir_path, artifact_name)

  # unpack artifact - creates the artifact directory
  artifact <- unpack_qiime2_artifact(qza, tmp_path)

  data_dir <- paste0(
    tmp_path, .Platform$file.sep,
    artifact$uuid, .Platform$file.sep,
    "data"
  )

  if (artifact$format == "TSVTaxonomyDirectoryFormat") {
    # bin taxonomy
    if (artifact$type == "FeatureData[Taxonomy]") {
      # read taxonomy table
      artifact[["data"]] <- read.table(file.path(data_dir, "taxonomy.tsv"),
        sep = "\t", header = TRUE, quote = ""
      )
      names(artifact[["data"]]) <- c("bin_names", "taxonomies", "confidence")

      if (remove_unpacked_artifacts) {
        unlink(tmp_path, recursive = TRUE)
      }

      return(artifact)
    }
  } else {
    if (remove_unpacked_artifacts) {
      unlink(tmp_path, recursive = TRUE)
    }
    abort_incorrect_format("TSVTaxonomyDirectoryFormat", artifact$format)
  }

  # only get here if it's a "TSVTaxonomyDirectoryFormat" but not
  # "FeatureData[Taxonomy]" type
  if (remove_unpacked_artifacts) {
    unlink(tmp_path, recursive = TRUE)
  }

  list()
}
