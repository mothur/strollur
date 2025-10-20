#' @title read_qiime2
#' @description
#' The read_qiime function reads various types of .qza files created by
#' \href{https://qiime2.org}{qiime2}, and creates a 'dataset' object.
#'
# nolint start
#' To generate the various input files you can follow \href{https://amplicon-docs.qiime2.org/en/latest/tutorials/moving-pictures.html}{qiime moving-pictures}.
# nolint end
#'
#' @param qza vector of filenames, .qza files containing your data from qiime2.
#' @param metadata filename, a .tsv file containing metadata
#' @param dataset_name A string containing a name for your dataset.
#' @param dir_path a string containing the name of directory where the artifacts
#' files should be unpacked. Default = current working directory.
#' @param remove_unpacked_artifacts boolean, When TRUE, the unpacked artifacts
#' and temporary directories will be removed. Default = TRUE.
#'
#' @examples
#'
#' # Using the example files from moving-pictures, we can assign bins, assign
#' # bin taxonomy and representative bin sequences, include a sample tree and
#' # metadata.
#'
#' qza_files <- c(
#'   rdataset_example("rep_seqs.qza"),
#'   rdataset_example("table.qza"),
#'   rdataset_example("taxonomy.qza"),
#'   rdataset_example("rooted-tree.qza")
#' )
#'
#' # data <- read_qiime2(
#' #   qza = qza_files,
#' #   metadata = rdataset_example("sample_metadata.tsv"),
#' #   dataset_name = "qiime_moving_pictures"
#' # )
#'
#' @return A 'dataset' object
#' @export
read_qiime2 <- function(qza, metadata = NULL,
                        dataset_name = "", dir_path = NULL,
                        remove_unpacked_artifacts = TRUE) {
  # if no dir_path given, set to current working directory
  if (is.null(dir_path)) {
    dir_path <- getwd()
  }

  if (!dir.exists(dir_path)) {
    # create it
    dir.create(dir_path)
  }

  data_found <- list()

  # extract data from the qza
  for (qza_file in qza) {
    artifact_name <- sub("\\.[^.]+$", "", basename(qza_file))
    tmp_path <- file.path(dir_path, artifact_name)

    # unpack artifact - creates the artifact directory
    artifact <- unpack_qiime2_artifact(qza_file, tmp_path)

    data_dir <- paste0(
      tmp_path, .Platform$file.sep,
      artifact$uuid, .Platform$file.sep,
      "data"
    )

    if (grepl("BIOMV", artifact$format)) {
      # feature table in biom format - shared bin data
      if (artifact$type == "FeatureTable[Frequency]") {
        # read biom file
        hdata <- read_biom(file.path(data_dir, "feature-table.biom"))

        # create data.frame from sparse otu data
        data_found[["bin_shared_assignments"]] <- data.frame(
          bin_names = hdata$otus[hdata$counts$i],
          abundances = hdata$counts$v,
          samples = hdata$samples[hdata$counts$j]
        )
      }
    } else if (
      (artifact$format == "DNASequencesDirectoryFormat") ||
        (artifact$format == "AlignedDNASequencesDirectoryFormat")
    ) {
      fasta_file <- "dna-sequences.fasta"
      if (artifact$format == "AlignedDNASequencesDirectoryFormat") {
        fasta_file <- "aligned-dna-sequences.fasta"
      }

      # bin representative sequences
      if (
        (artifact$type == "FeatureData[Sequence]") ||
          (artifact$type == "FeatureData[AlignedSequence]")
      ) {
        # read fasta file
        df <- read_fasta(file.path(data_dir, fasta_file))
        names(df) <- c("bin_names", "sequences")

        data_found[["bin_representatives"]] <- df
      }
    } else if (artifact$format == "TSVTaxonomyDirectoryFormat") {
      # bin taxonomy
      if (artifact$type == "FeatureData[Taxonomy]") {
        # read taxonomy table
        df <- read.table(file.path(data_dir, "taxonomy.tsv"),
          sep = "\t", header = TRUE, quote = ""
        )
        df <- df[, -c(3)]
        names(df) <- c("bin_names", "taxonomies")
        data_found[["bin_taxonomy"]] <- df
      }
    } else if (artifact$format == "NewickDirectoryFormat") {
      data_found[["sample_tree"]] <- ape::read.tree(file.path(
        data_dir,
        "tree.nwk"
      ))
    }

    if (remove_unpacked_artifacts) {
      unlink(tmp_path, recursive = TRUE)
    }
  }

  # read metadata
  if (!is.null(metadata)) {
      data_found[["metadata"]] <- read_qiime2_metadata(metadata)
  }

  # create new blank dataset
  data <- dataset$new(name = dataset_name)

  # add data_found to dataset in order that does not cause errors
  if (length(data_found) != 0) {

      data_names <- names(data_found)

      # if data includes shared data and representative sequences, assign
      # sample frequency to representative sequence and assign representative
      # sequence to bin
      if (all(c("bin_shared_assignments",
                "bin_representatives") %in% data_names) &&
          !(c("sequence_data") %in% data_names)) {

          data$add_sequences(data_found[["bin_representatives"]],
                                         sequence_names = "bin_names")
          data$assign_sequence_abundance(data_found[["bin_shared_assignments"]],
                                         sequence_names = "bin_names")
          data$assign_bins(bin_names = data_found[["bin_shared_assignments"]][,1],
                           sequence_names = data_found[["bin_shared_assignments"]][,1])
          data$assign_bin_representative_sequences(bin_names = data_found[["bin_shared_assignments"]][,1],
                           sequence_names = data_found[["bin_shared_assignments"]][,1])

      }

      # add sequence data

      # add sequence abundance

      # add sequence taxonomy

      # add sequence tree

      # add bin data

      # add bin taxonomy
      if ("bin_taxonomy" %in% data_names) {
          data$assign_bin_taxonomy(data_found[["bin_taxonomy"]])
      }

      # add sample tree
      if ("sample_tree" %in% data_names) {
          data$add_sample_tree(data_found[["sample_tree"]])
      }

      # add metadata

  }

  data
}
