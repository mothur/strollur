#' @title read_qiime2_feature_table
#' @description
#' Read a \href{https://qiime2.org}{qiime2} qza containing bin data
#'
#' # nolint start
#' To generate the various input files you can follow \href{https://amplicon-docs.qiime2.org/en/latest/tutorials/moving-pictures.html}{qiime moving-pictures}.
#' # nolint end
#'
#' @param qza file name, a qiime2 .qza file containing bin data.
#' @param dir_path a string containing the name of directory where the artifacts
#' files should be unpacked. Default = current working directory.
#' @param remove_unpacked_artifacts boolean, When TRUE, the artifact's temporary
#' directories will be removed after processing. Default = TRUE.
#'
#' @examples
#'
#' artifact <- read_qiime2_feature_table(rdataset_example(
#'   "table.qza"
#' ))
#'
#' # access the bin assignment table
#'
#' artifact$data
#'
#' # to create a 'dataset' object with your data
#'
#' data <- dataset$new("my_data")
#' data$assign_bins(artifact$data)
#' data
#'
#' @return A list containing artifact
#' @export
read_qiime2_feature_table <- function(qza, dir_path = NULL,
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

    if (grepl("BIOMV", artifact$format)) {
        # bin assignments
        if (artifact$type == "FeatureTable[Frequency]") {
            # read biom file
            hdata <- read_biom(file.path(data_dir, "feature-table.biom"))

            # create data.frame from sparse otu data
            artifact[["data"]] <- data.frame(
                bin_names = hdata$otus[hdata$counts$i],
                abundances = hdata$counts$v,
                samples = hdata$samples[hdata$counts$j]
            )

            if (remove_unpacked_artifacts) {
                unlink(tmp_path, recursive = TRUE)
            }

            return(artifact)
        }
    } else {
        if (remove_unpacked_artifacts) {
            unlink(tmp_path, recursive = TRUE)
        }
        abort_incorrect_format("BIOMV210DirFmt", artifact$format)
    }

    # only get here if it's a "BIOMV210DirFmt" but not
    # "FeatureTable[Frequency]" type
    if (remove_unpacked_artifacts) {
        unlink(tmp_path, recursive = TRUE)
    }

    list()
}
