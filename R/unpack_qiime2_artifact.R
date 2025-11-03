#' @title unpack_qiime2_artifact
#' @description
#' The unpack_qiime2_artifact function reads .qza files created by
#' \href{https://qiime2.org}{qiime2}, and returns the artifact.
#'
# nolint start
#' To generate the various input files you can follow \href{https://amplicon-docs.qiime2.org/en/latest/tutorials/moving-pictures.html}{qiime moving-pictures}.
# nolint end
#'
#' @param qza filename, a .qza file containing artifact
#' @param dir_path a string containing the name of directory where the artifacts
#' files should be written. Default = current working directory.
#'
#' @examples
#'
#' # Using the example files from moving-pictures
#'
#' artifact <- unpack_qiime2_artifact(
#'   qza = rdataset_example("table.qza"),
#' )
#'
#' @return A unpacked qza artifact
#' @export
unpack_qiime2_artifact <- function(qza, dir_path = NULL) {
  # error checks
  if (!file.exists(qza)) {
    .abort_nonexistant_file(qza)
  }

  # if no dir given, set to current working directory
  if (is.null(dir_path)) {
    artifact_name <- sub("\\.[^.]+$", "", basename(qza))
    dir_path <- file.path(getwd(), artifact_name)
  }

  if (!dir.exists(dir_path)) {
    # create it
    dir.create(dir_path)
  }

  # decompress qza file
  unzip(qza, exdir = dir_path)

  # get list of files in artifact
  file_list <- unzip(qza, exdir = dir_path, list = TRUE)

  # we need the uuid to find the metadata.yaml file
  # 6a560288-898e-4c1d-92ac-dd8d7822dcc9/metadata.yaml
  artifact_uuid <- gsub("/..+", "", file_list$Name[1])
  metadata <- paste0(
    dir_path, .Platform$file.sep,
    artifact_uuid, .Platform$file.sep,
    "metadata.yaml"
  )

  # reads uuid, type and format of artifact
  artifact <- read_yaml(metadata)
  artifact$contents <- data.frame(files = file_list)
  version_file <- paste0(
    dir_path, .Platform$file.sep,
    artifact_uuid, .Platform$file.sep,
    "VERSION"
  )
  artifact$version <- readr::read_table(version_file, show_col_types = FALSE)

  prov_names <- grep("..+provenance/..+action.yaml",
    file_list$Name,
    value = TRUE
  )
  # add provenance
  artifact$provenance <- lapply(paste0(
    dir_path,
    .Platform$file.sep,
    grep(paste0("..+provenance", .Platform$file.sep, "..+action.yaml"),
      file_list$Name,
      value = TRUE
    )
  ), read_yaml)

  names(artifact$provenance) <- prov_names

  artifact
}
