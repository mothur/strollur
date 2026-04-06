#' @title The `strollur` object stores the data associated with your microbial
#'   analysis.
#' @name strollur
#' @rdname strollur
#' @description 'strollur' is an R6 class that stores nucleotide sequences,
#' abundance, sample and treatment assignments, taxonomic classifications,
#' asv / otu clusters and various reports. It is designed to facilitate data
#' analysis across multiple R packages.
#'
#' @author Sarah Westcott, \email{swestcot@@umich.edu}
#'
#' @importFrom R6 R6Class
#' @importFrom methods new
#' @importFrom parallelly availableCores
#' @importFrom waldo compare
#' @import cli
#' @export
strollur <- R6Class("strollur",
  public = list(
    #' @field data Rcpp::XPtr<Dataset> pointer to 'Dataset' c++ class. This
    #' allows package developers an easy access point to the underlying C++ code
    #' with additional functionality.
    data = NULL,

    #' @field raw Rcpp::RawVector containing the serialized data of the
    #' 'Dataset' c++ class. This allows the load and save functions to work with
    #' the class.
    raw = NULL,

    #' @field sequence_tree a tree that relates sequences to eachother
    sequence_tree = NULL,

    #' @field sample_tree a tree that relates samples to eachother
    sample_tree = NULL,

    #' @description
    #' Create a new strollur dataset
    #' @param name String, name of dataset (optional)
    #' @param processors Integer, number of cores to use.
    #'  Default = all available
    #' @param dataset a `strollur` object.
    #' @examples
    #'
    #' # to create an empty strollur object, run the following:
    #'
    #' data <- new_dataset("soil")
    #'
    #' @return A new `strollur` object.
    initialize = function(name = "",
                          processors = parallelly::availableCores(),
                          dataset = NULL) {
      if (is.null(dataset)) {
        self$data <- xint_new_pointer(name, processors)
        self$sequence_tree <- NULL
        private$processors <- processors
        self$sample_tree <- NULL
      } else {
        # copy of dataset backend
        self$data <- xint_copy_pointer(dataset)
        xdev_set_num_processors(self, processors)
        # assign new name
        if (name != "") {
          xdev_set_dataset_name(self, name)
        }

        private$processors <- processors
        self$sequence_tree <- dataset$get_sequence_tree()
        self$sample_tree <- dataset$get_sample_tree()
      }

      invisible(self)
    },

    #' @description
    #' Get summary of `strollur` object
    print = function() {
      if (names(self, type = "dataset")[1] != "") {
        cat(names(self, type = "dataset")[1])
        cat(":\n\n")
      }

      # get dataset summaries
      results <- private$get_summary()

      results_names <- names(results)

      # print sequence summary - converting to ints
      if ("sequence_summary" %in% results_names) {
        for (x in 1:6) {
          results[["sequence_summary"]][, x] <- as.integer(
            results[["sequence_summary"]][, x]
          )
        }
        results[["sequence_summary"]]$numseqs <- sprintf(
          "%.2f", results[["sequence_summary"]]$numseqs
        )
        print(results[["sequence_summary"]])
      }

      if ("scrap_summary" %in% results_names) {
        cat("scrap_summary:\n")
        print(results[["scrap_summary"]])
      }

      if (count(data = self, type = "sequences", distinct = TRUE) != 0) {
        cat(
          paste("\nNumber of unique seqs:", count(
            data = self,
            type = "sequences",
            distinct = TRUE
          )),
          "\n"
        )
      } else {
        cat("\n")
      }
      cat(
        paste("Total number of seqs:", count(data = self, type = "sequences")),
        "\n\n"
      )

      # print number of samples
      if (count(data = self, type = "samples") != 0) {
        cat(paste0(
          "Total number of samples: ",
          count(data = self, type = "samples")
        ), "\n")
      }
      # print number of treatments
      if (count(data = self, type = "treatments") != 0) {
        cat(paste0(
          "Total number of treatments: ",
          count(data = self, type = "treatments")
        ), "\n")
      }

      # print number of each bin type
      bin_types <- get_bin_types(self)
      for (bin_type in bin_types) {
        if (count(data = self, type = "bins", bin_type = bin_type) != 0) {
          cat(
            paste0(
              "Total number of ", bin_type, "s: ",
              count(data = self, type = "bins", bin_type = bin_type)
            ),
            "\n"
          )
        }
      }
      # print number of resource references
      if (count(data = self, type = "references") != 0) {
        cat(paste0(
          "Total number of resource references: ",
          count(data = self, type = "references")
        ), "\n")
      }

      exclude <- c("sequence_scrap", "bin_scrap")
      report_names <- names(data = self, type = "reports")
      custom_report_names <- report_names[!report_names %in% exclude]

      if (length(custom_report_names) != 0) {
        cat(paste0(
          "Total number of custom reports: ",
          length(custom_report_names)
        ), "\n")
      }
      if (nrow(report(self, type = "metadata")) != 0) {
        cat(paste0("Your dataset includes metadata"), "\n")
      }
      cat("\n")
    },

    #' @description
    #' Get the abundance data for sequences, bins, samples, and treatments.
    #'
    #' @param type, string containing the type of data you want the number of.
    #'   Options include: "sequences", "bins", "samples" and "treatments".
    #'   Default = "sequences".
    #'
    #' @param bin_type, string containing the bin type you would like the
    #'   abundance data for. Default = "otu".
    #'
    #' @param by_sample, Boolean. When by_sample is TRUE, the abundance data
    #'   will be parsed by sample. Default = FALSE.
    #'
    #' @examples
    #'
    #' miseq <- miseq_sop_example()
    #'
    #' # To the total abundance for each sequence
    #' miseq$abundance(type = "sequences") |> head(n = 5)
    #'
    #' # To the total abundance for each sequence parsed by sample
    #' miseq$abundance(type = "sequences", by_sample = TRUE) |> head(n = 5)
    #'
    #' # To the total abundance for each "otu" bin
    #' miseq$abundance(type = "bins", bin_type = "otu") |> head(n = 5)
    #'
    #' # To the total abundance for each "otu" bin parsed by sample
    #' miseq$abundance(type = "bins", bin_type = "otu", by_sample = TRUE) |>
    #' head(n = 5)
    #'
    #' # To the total abundance for each "asv" bin
    #' miseq$abundance(type = "bins", bin_type = "asv") |> head(n = 5)
    #'
    #' # To the total abundance for each "asv" bin parsed by sample
    #' miseq$abundance(type = "bins", bin_type = "asv", by_sample = TRUE) |>
    #' head(n = 5)
    #'
    #' # To the total abundance for each sample
    #' miseq$abundance(type = "samples") |> head(n = 5)
    #'
    #' # To the total abundance for each treatment
    #' miseq$abundance(type = "treatments")
    abundance = function(type = "sequences",
                         bin_type = "otu",
                         by_sample = FALSE) {
      xdev_abundance(self, type, bin_type, by_sample)
    },

    #' @description
    #' Add sequences, reports, metadata or resource references
    #'
    #' @param table, a data.frame containing the data you wish to add.
    #'
    #' @param type, a string containing the type of data. Options include:
    #' 'sequences', 'references' 'metadata' and 'reports'.
    #'
    #' @param report_type, a string containing the type of report you are
    #' adding. Options include: 'metadata' and custom reports.
    #'
    #' @param table_names, named list used to indicate the names of the columns
    #'  in the table. By default:
    #'
    #' table_names <- list(sequence_name = "sequence_names",
    #'                     comment = "comments",
    #'                     sequence = "sequences",
    #'                     reference_name = "reference_names",
    #'                     reference_version = "reference_versions",
    #'                     reference_usage = "reference_usages",
    #'                     reference_note = "reference_notes",
    #'                     reference_url = "reference_urls")
    #'
    #' In table_names, 'sequence_name' is a string containing the name of the
    #' column in 'table' that contains the sequence names. It is used when you
    #' are adding FASTA data. Default column name is 'sequence_names'.
    #'
    #' In table_names, 'sequence' is a string containing the name of the column
    #' in 'table' that contains the sequence nucleotide strings. It is used when
    #' you are adding FASTA data. Default column name is 'sequences'.
    #'
    #' In table_names, 'comment' is a string containing the name of the column
    #' in 'table' that contains the sequence comments. It is used when you are
    #' adding FASTA data. Default column name is 'comments'.
    #'
    #' In table_names, 'reference_name' is a string containing the name of the
    #' column in 'table' that contains the reference names. It is used when you
    #' are adding reference data. Default column name is 'reference_names'.
    #'
    #' In table_names, 'reference_version' is a string containing the name of
    #' the column in 'table' that contains the reference versions. Default
    #' column name is 'reference_versions'.
    #'
    #' In table_names, 'reference_usage' is a string containing the name of the
    #' column in 'table' that contains the reference usages. Default column name
    #' is 'reference_usages'.
    #'
    #' In table_names, 'reference_note' is a string containing the name of the
    #' column in 'table' that contains the reference notes. Default column name
    #' is 'reference_notes'.
    #'
    #' In table_names, 'reference_url' is a string containing the name of the
    #' column in 'table' that contains the reference urls. Default column name
    #' is 'reference_urls'.
    #'
    #' @param reference, a list created by the function [new_reference].
    #'   Optional.
    #'
    #' @param verbose, boolean indicating whether or not you want progress
    #'   messages. Default = TRUE.
    #'
    #' @examples
    #'
    #' fasta_data <- read_fasta(fasta = strollur_example("final.fasta.gz"))
    #' contigs_report <- readRDS(strollur_example("miseq_contigs_report.rds"))
    #'
    #' # Create a new empty `strollur` object named 'example_dataset'
    #' data <- new_dataset(dataset_name = "example_dataset")
    #'
    #' data$add(table = fasta_data, type = "sequences")
    #' data$add(
    #'   table = contigs_report, type = "reports",
    #'   report_type = "contigs_report", list(sequence_name = "Name")
    #' )
    #'
    #' # To add metadata related to your study
    #'
    #' metadata <- readRDS(strollur_example("miseq_metadata.rds"))
    #'
    #' data$add(table = metadata, type = "metadata")
    #'
    #' # To add FASTA data with a resource reference
    #'
    #' # Create a new empty strollur object named 'example_dataset'
    #' data <- new_dataset(dataset_name = "example_dataset")
    #'
    #' # Create a resource reference for the FASTA data
    #' resource_url <- "https://mothur.org/wiki/silva_reference_files/"
    #' resource_reference <-
    #'   new_reference(
    #'     reference_name = "silva.bacteria.fasta",
    #'     reference_version = "1.38.1",
    #'     reference_usage = "alignment by mothur2 v1.0",
    #'     reference_note = "default options",
    #'     reference_url = resource_url
    #'   )
    #'
    #' # Add FASTA data with a resource reference
    #'
    #' data$add(
    #'   table = fasta_data,
    #'   type = "sequences",
    #'   reference = resource_reference
    #' )
    #'
    add = function(table,
                   type = "sequences",
                   report_type = NULL,
                   table_names = list(
                     sequence_name = "sequence_names",
                     sequence = "sequences",
                     comment = "comments",
                     reference_name = "reference_names",
                     reference_version = "reference_versions",
                     reference_usage = "reference_usages",
                     reference_note = "reference_notes",
                     reference_url = "reference_urls"
                   ),
                   reference = NULL,
                   verbose = TRUE) {
      add(self,
        table = table, type = type,
        report_type = report_type, table_names = table_names,
        reference = reference, verbose = verbose
      )

      invisible(self)
    },

    #' @description
    #' Add phylo tree relating the samples in your dataset
    #'
    #' @param tree a phylo tree object created by ape::read.tree.
    #' @examples
    #'
    #'  data <- strollur$new("my_dataset")
    #'
    #'  df <- read_mothur_shared(strollur_example("final.opti_mcc.shared"))
    #'  assign(data = data, table = df, type = "bins", bin_type = "otu")
    #'
    #'  tree <- ape::read.tree(strollur_example(
    #'  "final.opti_mcc.jclass.ave.tre"))
    #'
    #'  data$add_sample_tree(tree)
    #'
    add_sample_tree = function(tree) {
      if (!inherits(tree, "phylo")) {
        .abort_incorrect_type("phylo", tree)
      }

      # if no samples, add sequences in tree to dataset
      if (count(self, type = "samples") == 0) {
        message <- paste0("[Warning]: Your dataset does not contain sample ",
          "data, ignoring sample tree.",
          collapse = ""
        )
        cli::cli_alert(message)
      } else {
        # make sure the tree includes all "good" samples
        if (identical(
          sort(tree$tip.label),
          sort(names(self, type = "samples"))
        )) {
          # save tree
          self$sample_tree <- tree
        } else {
          # samples in dataset and not in tree
          missing_samples <- setdiff(
            names(self, type = "samples"),
            tree$tip.label
          )

          # if tree is "missing" names, then ignore tree
          if (length(missing_samples) != 0) {
            message <- paste("[WARNING]: Your tree does not",
              "contain a node for every sample in",
              "your dataset, ignoring tree.",
              "Missing tree nodes for:",
              paste(missing_samples, collapse = ", "),
              ".",
              collapse = ""
            )
            cli_alert(message)
          } else {
            # samples in tree and not in dataset
            extra_samples <- setdiff(
              tree$tip.label,
              names(self, type = "samples")
            )

            # if tree contains "extra" names, prune the tree
            self$sample_tree <- drop.tip(tree, tip = extra_samples)
          }
        }
      }

      invisible(self)
    },

    #' @description
    #' Add phylo tree relating the sequences in your dataset
    #' @param tree a phylo tree object created by ape::read.tree.
    #' @examples
    #'
    #'  data <- strollur$new("my_dataset")
    #'  tree <- ape::read.tree(strollur_example("final.phylip.tre.gz"))
    #'  data$add_sequence_tree(tree)
    #'
    add_sequence_tree = function(tree) {
      if (!inherits(tree, "phylo")) {
        .abort_incorrect_type("phylo", tree)
      }

      # if no seqs yet, add sequences in tree to dataset
      if (count(self, type = "sequences") == 0) {
        xdev_add_sequences(self, data.frame(sequence_names = tree$tip.label))

        # save tree
        self$sequence_tree <- tree
      } else {
        # make sure the tree includes all "good" sequences
        if (identical(
          sort(tree$tip.label),
          sort(names(self, type = "sequences"))
        )) {
          # save tree
          self$sequence_tree <- tree
        } else {
          # seqs in dataset and not in tree
          missing_seqs <- setdiff(
            names(self, type = "sequences"),
            tree$tip.label
          )

          # if tree is "missing" names, then ignore tree
          if (length(missing_seqs) != 0) {
            message <- paste("[WARNING]: Your tree does not",
              "contain a node for every sequence in",
              "your dataset, ignoring tree.",
              "Missing tree nodes for:",
              paste(missing_seqs, collapse = ", "),
              ".",
              collapse = ""
            )
            cli_alert(message)
          } else {
            # seqs in tree and not in dataset
            extra_seqs <- setdiff(
              tree$tip.label,
              names(self, type = "sequences")
            )

            # if tree contains "extra" names, prune the tree
            self$sequence_tree <- drop.tip(tree, tip = extra_seqs)
          }
        }
      }

      invisible(self)
    },

    #' @description
    #' Assign sequence abundances, sequence classifications, bins, bin
    #' representative sequences, bin classifications or treatments.
    #'
    #' @param table, a data.frame containing the data you wish to assign
    #'
    #' @param type, a string containing the type of data. Options include:
    #' 'sequence_abundance', 'sequence_taxonomy', 'bins',
    #'  'bin_representatives', 'bin_taxonomy' and 'treatments'.
    #'  Default = "bins".
    #'
    #' @param bin_type, string containing the bin type you would like the number
    #'   of bins for. Default = "otu".
    #'
    #' @param table_names, named list used to indicate the names of the columns
    #'   in the table. By default:
    #'
    #'   table_names <- list(sequence_name = "sequence_names", abundance =
    #'   "abundances", sample = "samples", treatment = "treatments", taxonomy =
    #'   "taxonomies", bin_name = "bin_names")
    #'
    #'   In table_names, 'sequence_name' is a string containing the name of the
    #'   column in 'table' that contains the sequence names. Default column name
    #'   is 'sequence_names'.
    #'
    #'   In table_names, 'abundance' is a string containing the name of the
    #'   column in 'table' that contains the abundances. Default column name is
    #'   'abundances'.
    #'
    #'   In table_names, 'sample' is a string containing the name of the column
    #'   in 'table' that contains the samples. Default column name is 'samples'.
    #'
    #'   In table_names, 'treatment' is a string containing the name of the
    #'   column in 'table' that contains the treatment names. Default column
    #'   name is 'treatments'.
    #'
    #'   In table_names, 'taxonomy' is a string containing the name of the
    #'   column in 'table' that contains the classifications. Default column
    #'   name is 'taxonomies'.
    #'
    #'   In table_names, 'bin_name' is a string containing the name of the
    #'   column in 'table' that contains the bin names. Default column name is
    #'   'bin_names'.
    #'
    #' @param reference, a list created by the function [new_reference].
    #'   Optional.
    #' @param verbose, boolean indicating whether or not you want progress
    #'   messages. Default = TRUE.
    #'
    #' @examples
    #'
    #' # create a new empty strollur object named 'example_dataset'
    #'
    #' data <- new_dataset(dataset_name = "example_dataset")
    #'
    #' # Assign sequence abundances
    #'
    #' abundance_by_sample <- read_mothur_count(strollur_example(
    #'   "final.count_table.gz"
    #' ))
    #'
    #' data$assign(table = abundance_by_sample, type = "sequence_abundance")
    #'
    #' # Assign sequence classifications
    #'
    #' sequence_classifications <- read_mothur_taxonomy(strollur_example(
    #'   "final.taxonomy.gz"
    #' ))
    #'
    #' data$assign(table = sequence_classifications, type = "sequence_taxonomy")
    #'
    #' # Assigning bins
    #'
    #' # read mothur's otu list file into data.frame
    #' otu_data <- read_mothur_list(list = strollur_example(
    #'   "final.opti_mcc.list.gz"
    #' ))
    #'
    #' # read mothur's asv list file into data.frame
    #' asv_data <- read_mothur_list(list = strollur_example(
    #'   "final.asv.list.gz"
    #' ))
    #'
    #' # read mothur's phylotype list file into data.frame
    #' phylo_data <- read_mothur_list(list = strollur_example(
    #'   "final.tx.list.gz"
    #' ))
    #'
    #' # read otu bin representative sequences into a data.frame
    #' bin_reps <- readRDS(strollur_example(
    #'                         "miseq_representative_sequences.rds"))
    #'
    #' # assign 'otu' bins using sequence names
    #' data$assign(table = otu_data, bin_type = "otu")
    #'
    #' # assign 'asv' bins using sequence names
    #' data$assign(table = asv_data, bin_type = "asv")
    #'
    #' # assign 'phylotype' bins using sequence names
    #' data$assign(table = phylo_data, bin_type = "phylotype")
    #'
    #' # assign 'otu' bin representative sequences
    #' data$assign(table = bin_reps, type = "bin_representatives")
    #'
    #' # To assign abundance only bins
    #'
    #' # create a new empty strollur object named 'example_dataset'
    #' data <- new_dataset(dataset_name = "example_dataset")
    #'
    #' # read mothur's shared file
    #' otu_data <- read_mothur_shared(strollur_example("final.opti_mcc.shared"))
    #'
    #' # assign abundance only otus parsed by sample
    #' data$assign(table = otu_data, bin_type = "otu")
    #'
    #' # Assigning bin classifications
    #'
    #' # read bin taxonomies
    #' otu_data <- read_mothur_cons_taxonomy(strollur_example(
    #'   "final.cons.taxonomy"
    #' ))
    #'
    #' # assign otu consensus taxonomies
    #' data$assign(
    #'   table = otu_data,
    #'   type = "bin_taxonomy", bin_type = "otu"
    #' )
    #'
    #' # Assign treatments
    #'
    #' sample_assignments <- readRDS(
    #'    strollur_example("miseq_sample_design.rds"))
    #'
    #' data$assign(table = sample_assignments, type = "treatments")
    #'
    #' @return double - The number of items assigned
    assign = function(table,
                      type = "bins",
                      bin_type = "otu",
                      table_names = list(
                        sequence_name = "sequence_names",
                        abundance = "abundances",
                        sample = "samples",
                        treatment = "treatments",
                        taxonomy = "taxonomies",
                        bin_name = "bin_names"
                      ),
                      reference = NULL,
                      verbose = TRUE) {
      assign(self,
        table = table, type = type,
        bin_type = bin_type, table_names = table_names,
        reference = reference, verbose = verbose
      )
      invisible(self)
    },

    #' @description
    #' Clear data from datasest
    clear = function() {
      clear(self)

      invisible(self)
    },

    #' @description
    #' Find the number of sequences, samples, treatments or bins of a given type
    #'
    #' @param type, string containing the type of data you want the number of.
    #' Options include: "sequences", "samples", "treatments", "bins", and
    #'  "references". Default = "sequences".
    #'
    #' @param bin_type, string containing the bin type you would like the number
    #'   of bins for. Default = "otu".
    #'
    #' @param samples, vector of strings. samples is only used when 'type' =
    #'   "sequences" or 'type' = "bins" . samples should contain the names of
    #'   the samples you want the count for. Default = NULL.
    #'
    #' @param distinct, Boolean. distinct is used when 'type' = "sequences" or
    #'   'type' = "bins". When 'type' = "sequences" and distinct is TRUE the
    #'   number of unique sequences is returned. When 'type' = "sequences" and
    #'   distinct is FALSE the total number of sequences is returned. This can
    #'   also be combined with samples to find the number of unique sequences
    #'   found ONLY in a given set of samples, or to find the number of unique
    #'   sequences in given set of samples that may also be present in other
    #'   samples. When 'type' = "bins", you can set distinct = TRUE to return
    #'   the number of bins that ONLY contain sequences from the given samples.
    #'   When distinct is FALSE the count returned contains bins with sequences
    #'   from a given samples, but those bins may also contain other samples.
    #'   Default = FALSE.
    #'
    #' @examples
    #'
    #' miseq <- miseq_sop_example()
    #'
    #' # To get the total number of sequences
    #' miseq$count(type = "sequences")
    #'
    #' # To get number of unique sequences
    #' miseq$count(type = "sequences", distinct = TRUE)
    #'
    #' # To get number of unique sequences from samples 'F3D0' and 'F3D1'
    #' # Note these sequences will be present in both samples but may be
    #' # be present in other samples as well
    #' miseq$count(type = "sequences", samples = c("F3D0", "F3D1"))
    #'
    #' # To get number of unique sequences exclusive to samples 'F3D0' and
    #' # 'F3D1'. Note sequences are present in both samples and NOT present in
    #' # any other samples.
    #'
    #' miseq$count(type = "sequences",
    #'             samples = c("F3D0", "F3D1"), distinct = TRUE )
    #'
    #' # To get the number of samples in the dataset
    #' miseq$count(type = "samples")
    #'
    #' # To get the number of treatments in the dataset
    #' miseq$count(type = "treatments")
    #'
    #' # To get the number of "otu" bins in the dataset
    #' miseq$count(type = "bins", bin_type = "otu")
    #'
    #' # To get the number of "asv" bins in the dataset
    #' miseq$count(type = "bins", bin_type = "asv")
    #'
    #' # To get the number of "phylotype" bins in the dataset
    #' miseq$count(type = "bins", bin_type = "phylotype")
    #'
    #' # To get number of "otu" bins from samples 'F3D0' and 'F3D1'
    #' # Note these bins will have sequences from both samples but there may be
    #' # other samples present as well
    #' miseq$count(
    #'   type = "bins", bin_type = "otu", samples = c("F3D0", "F3D1")
    #' )
    #'
    #' # To get number of "otu" bins unique to samples 'F3D0' and 'F3D1'
    #' # Note these bins will have sequences from both samples and NO other
    #' # samples will be present in the bins.
    #'
    #' miseq$count(
    #'   type = "bins", bin_type = "otu",
    #'   samples = c("F3D0", "F3D1"), distinct = TRUE
    #' )
    #'
    #' @return double
    count = function(type = "sequences",
                     bin_type = "otu",
                     samples = NULL,
                     distinct = FALSE) {
      xdev_count(self, type, bin_type, samples, distinct)
    },


    #' @description
    #' Get bin table types
    #' @examples
    #'
    #' data <- miseq_sop_example()
    #' data$get_bin_types()
    #'
    #' @return vector of strings
    get_bin_types = function() {
      get_bin_types(self)
    },

    #' @description
    #' Get phylo tree relating the samples in your dataset.
    #' @examples
    #'
    #'  tree <- ape::read.tree(strollur_example(
    #'   "final.opti_mcc.jclass.ave.tre"))
    #'
    #'  df <- read_mothur_shared(strollur_example("final.opti_mcc.shared"))
    #'
    #'  data <- strollur$new("my_dataset")
    #'
    #'  # assign abundance 'otu' bins
    #'  assign(data = data, table = df, type = "bins", bin_type = "otu")
    #'
    #'  data$add_sample_tree(tree)
    #'  data$get_sample_tree()
    #'
    get_sample_tree = function() {
      if (!is.null(self$sample_tree)) {
        # prune tree if needed
        # samples in tree and not in dataset
        extra_samples <- setdiff(
          self$sample_tree$tip.label,
          names(self, type = "samples")
        )

        if (length(extra_samples) != 0) {
          # if tree contains "extra" samples, prune the tree
          self$sample_tree <- drop.tip(self$sample_tree,
            tip = extra_samples
          )
        }
      }
      self$sample_tree
    },

    #' @description
    #' Get phylo tree relating the sequences in your strollur object.
    #' @examples
    #'
    #'  data <- strollur$new("my_dataset")
    #'  tree <- ape::read.tree(strollur_example("final.phylip.tre.gz"))
    #'  data$add_sequence_tree(tree)
    #'  data$get_sequence_tree()
    #'
    get_sequence_tree = function() {
      if (!is.null(self$sequence_tree)) {
        # prune tree if needed
        # seqs in tree and not in dataset
        extra_seqs <- setdiff(
          self$sequence_tree$tip.label,
          names(self, type = "sequences")
        )

        if (length(extra_seqs) != 0) {
          # if tree contains "extra" names, prune the tree
          self$sequence_tree <- drop.tip(self$sequence_tree,
            tip = extra_seqs
          )
        }
      }
      self$sequence_tree
    },

    #' @description
    #' Get the names of a given type of data
    #'
    #' @param type, string containing the type of data you would like. Options
    #'   include: "dataset", "sequences", "bins", "samples", "treatments",
    #'   "reports". Default = "sequences".
    #'
    #' @param bin_type, string containing the bin type you would like the names
    #'   for. Default = "otu".
    #'
    #' @param samples, vector of strings. samples is only used when 'type' =
    #'   "sequences" or 'type' = "bins" . samples should contain the names of
    #'   the samples you want names for. Default = NULL.
    #'
    #' @param distinct, Boolean. distinct is used when 'type' = "sequences" or
    #'   'type' = "bins" and the samples parameter is used. The distinct
    #'   parameter allows you to get the names that present given set of
    #'   samples. When distinct is TRUE, the names function will return the
    #'   names that ONLY contain data from the given samples. When distinct is
    #'   FALSE the data returned contains data from a given samples, but may
    #'   ALSO contain data from other samples. Default = FALSE.
    #'
    #' @examples
    #'
    #' miseq <- miseq_sop_example()
    #'
    #' # To get the name of the dataset
    #' miseq$names(type = "dataset")
    #'
    #' # To get the names of the sequences
    #' miseq$names(type = "sequences")
    #'
    #' # To get the names of the sequences present sample 'F3D0'
    #' miseq$names(type = "sequences", samples = c("F3D0"))
    #'
    #' #' # To get the names of the sequences unique to sample 'F3D0'
    #' miseq$names(type = "sequences", samples = c("F3D0"), distinct = TRUE)
    #'
    #' # To get the names of the samples
    #' miseq$names(type = "samples")
    #'
    #' # To get the names of the treatments
    #' miseq$names(type = "treatments")
    #'
    #' # To get the names of the bins
    #' miseq$names(type = "bins")
    #'
    #' # To get the names of the bins that are unique to 'F3D0'
    #' miseq$names(type = "bins", samples = c("F3D0"), distinct = TRUE)
    #'
    #' # To get the names of the bins that include sequences from 'F3D0'
    #' miseq$names(type = "bins", samples = c("F3D0"), distinct = FALSE)
    #'
    #' # To get the names of the reports
    #' miseq$names(type = "reports")
    #'
    #' @return vector of strings, containing the names requested
    names = function(type = "sequences",
                     bin_type = "otu",
                     samples = NULL,
                     distinct = FALSE) {
      xdev_names(self, type, bin_type, samples, distinct)
    },

    #' @description
    #' Get a data.frame containing the given report
    #'
    #' @param type, string containing the type of report you would like. Options
    #' include: "fasta", "sequences", "sequence_bin_assignments",
    #' "sequence_taxonomy", "bin_taxonomy", "bin_representatives",
    #'  "sample_assignments", "metadata", "references", "sequence_scrap",
    #' "bin_scrap". If you have added custom reports for alignment,
    #' contigs_assembly or chimeras, you can get those as well.
    #'  Default = "sequences".
    #'
    #' @param bin_type, string containing the bin type you would like a
    #'   bin_taxonomy report for. Default = "otu".
    #'
    #' @examples
    #'
    #' miseq <- miseq_sop_example()
    #'
    #' # To get the FASTA data
    #'
    #' miseq$report(type = "fasta") |> head(n = 5)
    #'
    #' # To get a report about the FASTA data
    #'
    #' miseq$report(type = "sequences") |> head(n = 5)
    #'
    #' # To get the sequence bin assignments
    #'
    #' miseq$report(type = "sequence_bin_assignments", bin_type = "otu") |>
    #' head(n = 5)
    #'
    #' # To get the sample treatment assignments
    #'
    #' miseq$report(type = "sample_assignments") |> head(n = 5)
    #'
    #' # To get a report about sequence classifications
    #'
    #' miseq$report(type = "sequence_taxonomy") |> head(n = 5)
    #'
    #' # To get a report about bin classifications for 'otu' data
    #'
    #' miseq$report(type = "bin_taxonomy", bin_type = "otu") |> head(n = 5)
    #'
    #' # To get the 'otu' bin representative sequences
    #'
    #' miseq$report(type = "bin_representatives", bin_type = "otu") |>
    #' head(n = 5)
    #'
    #' # To get a report about the sequences removed during your analysis:
    #'
    #' miseq$report(type = "sequence_scrap")
    #'
    #' # To get a report about the "otu" bins removed during your analysis:
    #'
    #' miseq$report(type = "bin_scrap", bin_type = "otu")
    #'
    #' # To get the metadata associated with your data:
    #'
    #' metadata <- miseq$report(type = "metadata") |> head(n = 5)
    #'
    #' # To get the resource references associated with your data:
    #'
    #' references <- miseq$report(type = "references")
    #'
    #' # To get our custom report containing the contigs assembly data:
    #'
    #' miseq$report(type = "contigs_report") |> head(n = 5)
    #'
    #' @return data.frame
    report = function(type = "sequences", bin_type = "otu") {
      xdev_report(self, type, bin_type)
    },

    #' @description
    #' Summarize the sequences data, custom reports, and scrapped data
    #'
    #' @param type, string containing the type of data you want the number of.
    #' Options include: "sequences", "reports" and "scrap".
    #'  Default = "sequences".
    #'
    #' @param report_type, string containing the report type you would
    #'  summarized. For example, the miseq_sop_example includes contigs assembly
    #'  data and can be accessed with report_type = "contigs_report".
    #'  Default = NULL.
    #'
    #' @param verbose, boolean indicating whether or not you want progress
    #'  messages. Default = TRUE.
    #'
    #' @examples
    #'
    #' miseq <- miseq_sop_example()
    #'
    #' # To get the summary of your FASTA data
    #' miseq$summary(type = "sequences")
    #'
    #' # summarize contigs_report
    #' miseq$summary(type = "reports", report_type = "contigs_report")
    #'
    #' # remove sample 'F3D0' to produce a scrap report
    #' xdev_remove_samples(data = miseq, samples = c("F3D0"))
    #'
    #' # summarize scrapped data -
    #' # sequences and bins scrapped by removing the sample "F3D0"
    #' miseq$summary(type = "scrap")
    #'
    #' @return data.frame
    summary = function(type = "sequences",
                       report_type = NULL, verbose = TRUE) {
      dataset_summary <- xdev_summarize(self, type, report_type)
      if (verbose) {
        print(dataset_summary)
      }
      dataset_summary
    }
  ),
  private = list(
    processors = 1,
    version = "1.0.0",
    finalize = function() {},


    #' Get summary of the sequence reports
    get_summary = function() {
      results <- list()

      # if you have summary results to print
      if (!all(xdev_get_sequences(self) == "")) {
        # if you have summary results to print
        results[["sequence_summary"]] <- xdev_summarize(
          data = self,
          type = "sequences"
        )
      }

      exclude <- c("sequence_scrap", "bin_scrap")
      report_names <- names(data = self, type = "reports")
      report_names <- report_names[!report_names %in% exclude]

      if (length(report_names) != 0) {
        for (name in report_names) {
          results[[name]] <- xdev_summarize(
            data = self,
            type = "reports",
            report_type = name
          )
        }
      }

      df <- xdev_summarize(
        data = self,
        type = "scrap"
      )
      if (nrow(df) != 0) {
        results[["scrap_summary"]] <- df
      }

      if (count(self, type = "samples") != 0) {
        results[["sample_summary"]] <- abundance(self, type = "samples")

        if (count(self, "treatments") != 0) {
          results[["treatment_summary"]] <- abundance(self, type = "treatments")
        }
      }

      results
    }
  )
)
