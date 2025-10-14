#' @title read_qiime2
#' @description
#' The read_qiime function reads various types of .qza files created by
#' \href{https://qiime2.org}{qiime2}, and creates a 'dataset' object.
#'
# nolint start
#' To generate the various input files you can follow \href{https://amplicon-docs.qiime2.org/en/latest/tutorials/moving-pictures.html}{qiime moving-pictures}.
# nolint end
#'
#' @param biom filename, a .qza file containing a feature table in biom
#' format. This will be used to assign bin abundance and sample data.
#' @param fasta filename, a .qza file containing DNA sequence strings.
#' @param cons_fasta filename, a .qza file containing DNA sequence strings of
#' representative sequences for each bin.
#' @param taxonomy filename, a .qza file containing classifications.
#' @param cons_taxonomy filename, a .qza file containing conscensus
#' classifications for each bin.
#' @param tree filename, a .qza file containing a tree.
#' @param metadata filename, a .tsv file containing metadata
#' @param dataset_name A string containing a name for your dataset.
#' @param dir_path a string containing the name of directory where the artifacts
#' files should be unpacked. Default = current working directory.
#' @param remove_unpacked_artifacts boolean, When TRUE, the unpacked artifacts
#' and temporary directories will be removed. Default = TRUE.
#'
#' @examples
#'
#' # Using the example files from moving-pictures
#'
#' # dataset <- read_qiime2(
#' #   cons_fasta = rdataset_example("rep_seqs.qza"),
#' #   biom = rdataset_example("table.qza"),
#' #   cons_taxonomy = rdataset_example("taxonomy.qza"),
#' #   tree = rdataset_example("rooted-tree.qza"),
#' #   metadata = rdataset_example("sample_metadata.tsv"),
#' #   dataset_name = "qiime_moving_pictures"
#' # )
#'
#' @return A 'dataset' object
#' @export
read_qiime2 <- function(fasta = NULL, taxonomy = NULL,
                        biom = NULL, cons_fasta = NULL, cons_taxonomy = NULL,
                        tree = NULL, metadata = NULL,
                        dataset_name = "", dir_path = NULL,
                        remove_unpacked_artifacts = TRUE) {
  # if no dir_path given, set to current working directory
  if (is.null(dir_path)) {
    dir_path <- getwd()
  }

  # create new blank dataset
  data <- dataset$new(name = dataset_name)

  artifact_directories_to_remove <- c()

  # add bin assignments - shared data
  if (!is.null(biom)) {
    artifact_name <- sub("\\.[^.]+$", "", basename(biom))
    tmp_path <- file.path(dir_path, artifact_name)

    # unpack artifact - creates the artifact directory
    artifact <- unpack_qiime2_artifact(biom, tmp_path)

    artifact_directories_to_remove <- c(
      artifact_directories_to_remove,
      tmp_path
    )

    feature_table_dir <- paste0(
      tmp_path, .Platform$file.sep,
      artifact$uuid, .Platform$file.sep,
      "data"
    )

    # read biom file, bin assignments by sample
    hdata <- read_biom(file.path(
      feature_table_dir,
      "feature-table.biom"
    ))

    # create data.frame from sparse otu data
    df <- data.frame(
      bin_names = hdata$otus[hdata$counts$i],
      abundances = hdata$counts$v,
      samples = hdata$samples[hdata$counts$j]
    )

    # add bin assignments to dataset
    data$assign_bins(df)
  }

  # add bin representative assignments
  if (!is.null(cons_fasta)) {
    artifact_name <- sub("\\.[^.]+$", "", basename(cons_fasta))
    tmp_path <- file.path(dir_path, artifact_name)

    # unpack artifact - creates the artifact directory
    artifact <- unpack_qiime2_artifact(cons_fasta, tmp_path)

    artifact_directories_to_remove <- c(
      artifact_directories_to_remove,
      tmp_path
    )

    # confirm these are DNA sequences
    if (
      (artifact$format != "DNASequencesDirectoryFormat") &&
        (artifact$format != "AlignedDNASequencesDirectoryFormat")
    ) {
      message <- paste0(
        "[ERROR]: ", fasta, " artifact format is ",
        artifact$format, ", expected formats are ",
        "'AlignedDNASequencesDirectoryFormat' and ",
        "'DNASequencesDirectoryFormat'."
      )
      cli::cli_abort(message)
    }

    fasta_dir <- paste0(
      tmp_path, .Platform$file.sep,
      artifact$uuid, .Platform$file.sep,
      "data"
    )

    fasta_file <- "dna-sequences.fasta"
    if (artifact$format == "AlignedDNASequencesDirectoryFormat") {
      fasta_file <- "aligned-dna-sequences.fasta"
    }

    # read fasta file
    df <- read_fasta(file.path(fasta_dir, fasta_file))
    names(df) <- c("bin_names", "sequences")

    # add sequences to dataset
    data$assign_bin_representative_sequences(df)
  }

  if (remove_unpacked_artifacts) {
    for (artifact_dir in artifact_directories_to_remove) {
      unlink(artifact_dir, recursive = TRUE)
    }
  }

  data
}
