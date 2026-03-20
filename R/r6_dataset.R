#' @title dataset
#' @description 'dataset' is an R6 class that stores nucleotide sequences,
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
dataset <- R6Class("dataset",
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
    #' Create a new dataset
    #' @param name String, name of dataset (optional)
    #' @param processors Integer, number of cores to use.
    #'  Default = all available
    #' @param dataset a `dataset` object.
    #' @examples
    #'
    #' # to create an empty dataset, run the following:
    #'
    #' data <- new_dataset("soil")
    #'
    #' @return A new `dataset` object.
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
    #' Get summary of sequence data
    print = function() {
      if (names(self, type = "dataset")[1] != "") {
        cat(names(self, type = "dataset")[1])
        cat(":\n\n")
      }

      # get dataset summaries
      results <- self$get_summary()

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
    #' Add phylo tree relating the samples in your dataset
    #'
    #' @param tree a phylo tree object created by ape::read.tree.
    #' @examples
    #'
    #'  data <- dataset$new("my_dataset")
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
    #'  data <- dataset$new("my_dataset")
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
    #' Clear data from datasest
    clear = function() {
      clear(self)

      invisible(self)
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
    #' Get data.frame containing metadata for the dataset
    #' @examples
    #'   data <- dataset$new("my_dataset")
    #'
    #'   metadata <- readr::read_tsv(strollur_example("sample-metadata.tsv"),
    #'    col_names = TRUE, show_col_types = FALSE)
    #'
    #'   add(data = data, table = metadata, type = "metadata")
    #'
    #'   data$get_metadata()
    #'
    #' @return data.frame()
    get_metadata = function() {
      report(self, "metadata")
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
    #'  data <- dataset$new("my_dataset")
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
    #' Get data.frame sequence report data. Sequence report data includes: start
    #' positions, end positions, number of bases, number of ambiguous bases,
    #' length of longest homopolymer, and the number of N's.
    #' @return data.frame
    get_sequence_report = function() {
      report(self, "sequences")
    },

    #' @description
    #' Get summary of the sequence reports
    #' @return list of data.frames
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
    },

    #' @description
    #' Get phylo tree relating the sequences in your dataset.
    #' @examples
    #'
    #'  data <- dataset$new("my_dataset")
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
    }
  ),
  private = list(
    processors = 1,
    version = "1.0.0",
    finalize = function() {}
  )
)
