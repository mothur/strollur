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
#' # Using the example files from moving-pictures, we add FASTA data, assign
#' # taxonomy and abundance for features, and add a newick tree and
#' # metadata.
#'
#' qza_files <- c(
#'   rdataset_example("rep_seqs.qza"),
#'   rdataset_example("table.qza"),
#'   rdataset_example("taxonomy.qza"),
#'   rdataset_example("rooted-tree.qza")
#' )
#'
#' data <- read_qiime2(
#'   qza = qza_files,
#'   metadata = rdataset_example("sample_metadata.tsv"),
#'   dataset_name = "qiime_moving_pictures"
#' )
#' data
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

  # types of data -> bin_shared_assignments, bin_representatives,
  #  bin_taxonomy, sequence_tree, metadata
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
      data_found[["sequence_tree"]] <- ape::read.tree(file.path(
        data_dir,
        "tree.nwk"
      ))
    }else {
        message <- paste0("Format ", artifact$format,
                          " is not supported, ignoring")
        cli_alert(message)
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
    bin_type <- "asv"

    # add fasta data
    if ("bin_representatives" %in% data_names) {
      add(data = data, table = data_found[["bin_representatives"]],
            type = "sequences", table_names = list(sequence_name = "bin_names"))
    }

    # assign sequence abundance by sample
    if ("bin_shared_assignments" %in% data_names) {
        assign(data = data, table = data_found[["bin_shared_assignments"]],
               type = "sequence_abundance",
               table_names = list(sequence_name = "bin_names"))
    }

    # assign sequences to bins
    if ("bin_shared_assignments" %in% data_names) {
        assign(data = data, table = data_found[["bin_shared_assignments"]],
               type = "bins", bin_type = bin_type)
    }

    # add bin taxonomy
    if ("bin_taxonomy" %in% data_names) {
      assign(data = data, table = data_found[["bin_taxonomy"]],
             type = "bin_taxonomy", bin_type = bin_type)
    }

    # add sequence tree
    if ("sequence_tree" %in% data_names) {
      data$add_sequence_tree(data_found[["sequence_tree"]])
    }
  }

  data
}
