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
#' @importFrom parallelly, availableCores
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
      self$get_summary()
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
        "\n"
      )

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
    #'  df <- read_mothur_shared(rdataset_example("final.opti_mcc.shared"))
    #'  assign_bins(data, df)
    #'
    #'  tree <- ape::read.tree(rdataset_example(
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
    #'  tree <- ape::read.tree(rdataset_example("final.phylip.tre"))
    #'  data$add_sequence_tree(tree)
    #'
    add_sequence_tree = function(tree) {
      if (!inherits(tree, "phylo")) {
        .abort_incorrect_type("phylo", tree)
      }

      # if no seqs yet, add sequences in tree to dataset
      if (count(self, type = "sequences") == 0) {
        add_sequences(self, data.frame(sequence_names = tree$tip.label))

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
    #' Remove 'sample_tree', or 'sequence_tree' data from your dataset.
    #' @param tags a vector of strings containing the items you wish to clear.
    #' Options are 'metadata', 'references', 'sequence_tree', 'sample_tree',
    #' 'alignment_report', 'contigs_assembly_report' and '"'chimera_report'.
    #' By default, everything is cleared.
    clear = function(tags = NULL) {
      if (is.null(tags)) {
        tags <- ""
      }

      if (tags == "") {
        self$sequence_tree <- NULL
        self$sample_tree <- NULL
      }

      valid_tags <- c("sequence_tree", "sample_tree")

      for (tag in tags) {
        if (!(tag %in% valid_tags)) {
          cli_alert("{.var {tag}} is not a valid item to clear, ignoring.")
        }
      }

      if ("sequence_tree" %in% tags) {
        self$sequence_tree <- NULL
      }
      if ("sample_tree" %in% tags) {
        self$sample_tree <- NULL
      }

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
    #'   metadata <- readr::read_tsv(rdataset_example("sample-metadata.tsv"),
    #'    col_names = TRUE, show_col_types = FALSE)
    #'
    #'   add_report(data, metadata, "metadata")
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
    #'  tree <- ape::read.tree(rdataset_example(
    #'   "final.opti_mcc.jclass.ave.tre"))
    #'
    #'  df <- read_mothur_shared(rdataset_example("final.opti_mcc.shared"))
    #'
    #'  data <- dataset$new("my_dataset")
    #'  assign_bins(data, df)
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
      report(self, "sequence_data")
    },

    #' @description
    #' Get summary of the sequence reports
    #' @param silent Default = FALSE, meaning print summaries
    #' @return list of data.frames
    get_summary = function(silent = FALSE) {
      results <- get_sequence_summary(self)

      non_report_names <- c("sequence_summary", "scrap_summary")
      t <- names(results)
      report_names <- t[!t %in% non_report_names]

      # if you have summary results to print
      if (!all(get_sequences(self) == "") && (!silent)) {
        # if you have summary results to print
        if (!silent) {
          cat("sequence_summary:\n")
          print(results[["sequence_summary"]])
          cat("Unique seqs:\t", count(
            data = self,
            type = "sequences", distinct = TRUE
          ), "\n")
          cat("Total seqs:\t", count(data = self, type = "sequences"), "\n")
        }
      }

      if (length(report_names) != 0) {
        if (!silent) {
          cat("\n")
          for (name in report_names) {
            cat(name, ":\n")
            print(results[[name]])
            cat("Unique seqs:\t", count(
              data = self,
              type = "sequences", distinct = TRUE
            ), "\n")
            cat("Total seqs:\t", count(data = self, type = "sequences"), "\n\n")
          }
        }
      }

      if ("scrap_summary" %in% names(results)) {
        if (!silent) {
          cat("\nTrash_code   Unique_count    Total_count:\n")
          for (i in seq_along(results$scrap_summary$trash_codes)) {
            cat(paste(
              results$scrap_summary$trash_codes[i],
              results$scrap_summary$unique_count[i],
              results$scrap_summary$total_count[i],
              sep = "\t"
            ), "\n")
          }
        }
      }

      if (count(self, type = "samples") != 0) {
        results[["sample_summary"]] <- totals(self, type = "samples")

        if (!silent) {
          cat("\nSample   Total:\n")
          sample_names <- results[["sample_summary"]]$samples
          sample_totals <- results[["sample_summary"]]$totals

          for (i in seq_along(sample_names)) {
            cat(paste(sample_names[i], sample_totals[i], sep = "\t"), "\n")
          }
        }

        if (count(self, "treatments") != 0) {
          results[["treatment_summary"]] <- totals(self, "treatments")

          if (!silent) {
            treatment_names <- results[["treatment_summary"]]$treatments
            treatment_totals <- results[["treatment_summary"]]$totals
            cat("\n")
            cat("Treatment   Total:\n")
            for (i in seq_along(treatment_names)) {
              cat(
                paste(treatment_names[i], treatment_totals[i],
                  sep = "\t"
                ), "\n"
              )
            }
          }
        }
      } else {
        if (!silent) {
          cli::cli_alert("Your dataset does not include sample data, ignoring.")
        }
      }

      return(results)
    },

    #' @description
    #' Get phylo tree relating the sequences in your dataset.
    #' @examples
    #'
    #'  data <- dataset$new("my_dataset")
    #'  tree <- ape::read.tree(rdataset_example("final.phylip.tre"))
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
