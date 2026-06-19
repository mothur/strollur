#' @title The `strollur` object stores the data associated with your amplicon
#'   sequence analysis.
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

    #' @field sequence_tree a tree that relates sequences to each other
    sequence_tree = NULL,

    #' @field sample_tree a tree that relates samples to each other
    sample_tree = NULL,

    #' @description
    #' Create a new strollur dataset
    #' @param name String, name of dataset (optional)
    #' @param dataset a `strollur` object.
    #' @examples
    #'
    #' # to create an empty strollur object, run the following:
    #'
    #' data <- new_dataset("soil")
    #'
    #' @return A new `strollur` object.
    initialize = function(name = "",
                          dataset = NULL) {
      if (is.null(dataset)) {
        self$data <- xint_new_pointer(dataset_name = name)
        self$sequence_tree <- NULL
        self$sample_tree <- NULL
      } else {
        if (!inherits(dataset, "strollur")) {
          stop("dataset must be a strollur object.")
        }

        # copy of dataset backend
        self$data <- xint_copy_pointer(dataset)
        # assign new name
        if (name != "") {
          xdev_set_dataset_name(self, name)
        }

        self$sequence_tree <- dataset$get_sequence_tree()
        self$sample_tree <- dataset$get_sample_tree()
      }

      invisible(self)
    },

    #' @description
    #' Print summary of `strollur` object
    #' @examples
    #' miseq <- load_dataset(strollur_example("miseq_sop.rds"))
    #' miseq
    #'
    #' @return No return value, called for side effects.
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
        cat("\nscrap_summary:\n")
        print(results[["scrap_summary"]])
      }

      if (xdev_count(data = self, type = "sequence", distinct = TRUE) != 0) {
        cat(
          paste("\nNumber of unique seqs:", count(
            data = self,
            type = "sequence",
            distinct = TRUE
          )),
          "\n"
        )
      } else {
        cat("\n")
      }
      cat(paste(
        "Total number of seqs:",
        xdev_count(data = self, type = "sequence")
      ), "\n\n")

      # print number of samples
      if (xdev_count(data = self, type = "sample") != 0) {
        cat(paste0(
          "Total number of samples: ",
          xdev_count(data = self, type = "sample")
        ), "\n")
      }
      # print number of treatments
      if (xdev_count(data = self, type = "treatment") != 0) {
        cat(paste0(
          "Total number of treatments: ",
          xdev_count(data = self, type = "treatment")
        ), "\n")
      }

      # print number of each bin type
      bin_types <- get_bin_types(self)
      for (bin_type in bin_types) {
        if (xdev_count(data = self, type = "bin", bin_type = bin_type) != 0) {
          cat(
            paste0(
              "Total number of ", bin_type, "s: ",
              xdev_count(data = self, type = "bin", bin_type = bin_type)
            ),
            "\n"
          )
        }
        bin_tax_report <- xdev_report(
          data = self,
          type = "bin_taxonomy",
          bin_type = bin_type
        )
        if (nrow(bin_tax_report) != 0) {
          cat(
            paste0(
              "Total number of ", bin_type, " bin classifications: ",
              length(unique(bin_tax_report[["bin_name"]]))
            ),
            "\n"
          )
        }
      }

      if (xdev_has_sequence_taxonomy(self)) {
        seq_tax_report <- xdev_report(
          data = self,
          type = "sequence_taxonomy"
        )
        cat(
          paste0(
            "Total number of sequence classifications: ",
            length(unique(seq_tax_report[["sequence_name"]]))
          ),
          "\n"
        )
      }

      # print number of resource references
      if (xdev_count(data = self, type = "resource_reference") != 0) {
        cat(paste0(
          "Total number of resource references: ",
          xdev_count(data = self, type = "resource_reference")
        ), "\n")
      }

      exclude <- c("sequence_scrap", "bin_scrap")
      report_names <- names(data = self, type = "report")
      custom_report_names <- report_names[!report_names %in% exclude]

      if (length(custom_report_names) != 0) {
        cat(paste0(
          "Total number of custom reports: ",
          length(custom_report_names)
        ), "\n")
      }
      cat("\n")
    },

    #' @description
    #' Get the abundance data for sequences, bins, samples, and treatments.
    #'
    #' @param type string containing the type of data you want the number of.
    #'   Options include: "sequence", "bin", "sample" and "treatment".
    #'   Default = "sequence".
    #'
    #' @param bin_type string containing the bin type you would like the
    #'   abundance data for. Default = "otu".
    #'
    #' @param by_sample Boolean. When by_sample is TRUE, the abundance data
    #'   will be parsed by sample. Default = FALSE.
    #'
    #' @examples
    #'
    #' miseq <- load_dataset(strollur_example("miseq_sop.rds"))
    #'
    #' # To the total abundance for each sequence
    #' miseq$abundance(type = "sequence") |> head(n = 5)
    #'
    #' # To the total abundance for each sequence parsed by sample
    #' miseq$abundance(type = "sequence", by_sample = TRUE) |> head(n = 5)
    #'
    #' # To the total abundance for each "otu" bin
    #' miseq$abundance(type = "bin", bin_type = "otu") |> head(n = 5)
    #'
    #' # To the total abundance for each "otu" bin parsed by sample
    #' miseq$abundance(type = "bin", bin_type = "otu", by_sample = TRUE) |>
    #' head(n = 5)
    #'
    #' # To the total abundance for each "asv" bin
    #' miseq$abundance(type = "bin", bin_type = "asv") |> head(n = 5)
    #'
    #' # To the total abundance for each "asv" bin parsed by sample
    #' miseq$abundance(type = "bin", bin_type = "asv", by_sample = TRUE) |>
    #' head(n = 5)
    #'
    #' # To the total abundance for each sample
    #' miseq$abundance(type = "sample") |> head(n = 5)
    #'
    #' # To the total abundance for each treatment
    #' miseq$abundance(type = "treatment")
    #'
    #' @return data.frame
    abundance = function(type = "sequence",
                         bin_type = "otu",
                         by_sample = FALSE) {
      xdev_abundance(self, type, bin_type, by_sample)
    },

    #' @description
    #' Add sequences, reports or resource references
    #'
    #' @param table a data.frame containing the data you wish to add.
    #' @param type a string containing the type of data. Options include:
    #' 'sequence', 'resource_reference' and 'report'.
    #' @param report_type a string containing the type of report you are
    #' adding. Options include: 'metadata' and custom reports.
    #' @param table_names named list used to indicate the names of the columns
    #'  in the table. By default:
    #'
    #' table_names <- list(sequence_name = "sequence_name",
    #'                     comment = "comment",
    #'                     sequence = "sequence",
    #'                     reference_name = "name",
    #'                     reference_vendor = "vendor",
    #'                     reference_version = "version",
    #'                     reference_usage = "usage",
    #'                     reference_note = "note",
    #'                     reference_documentation_url = "documentation_url",
    #'                     reference_method_url = "method_url",
    #'                     reference_parameter = "parameter",
    #'                     reference_citation = "citation")
    #'
    #' In table_names, 'sequence_name' is a string containing the name of the
    #' column in 'table' that contains the sequence names. It is used when you
    #' are adding FASTA data. Default column name is 'sequence_name'.
    #'
    #' In table_names, 'sequence' is a string containing the name of the column
    #' in 'table' that contains the sequence nucleotide strings. It is used when
    #' you are adding FASTA data. Default column name is 'sequence'.
    #'
    #' In table_names, 'comment' is a string containing the name of the column
    #' in 'table' that contains the sequence comments. It is used when you are
    #' adding FASTA data. Default column name is 'comment'.
    #'
    #' In table_names, 'reference_vendor' is a string containing the name of the
    #' column in 'table' that contains the reference vendor names. It is used
    #' when ' you are adding reference data. Default column name is 'vendor'.
    #'
    #' In table_names, 'reference_name' is a string containing the name of the '
    #' column in 'table' that contains the reference names. It is used when you
    #' are ' adding reference data. Default column name is 'name'.
    #'
    #' In table_names, 'reference_version' is a string containing the name of
    #' the ' column in 'table' that contains the reference versions. Default
    #' column name is 'version'.
    #'
    #' In table_names, 'reference_usage' is a string containing the name of the
    #' column in 'table' that contains the reference usages. Default column name
    #' is 'usage'.
    #'
    #' In table_names, 'reference_note' is a string containing the name of the
    #' column in 'table' that contains the reference notes. Default column name
    #' is 'note'.
    #'
    #' In table_names, 'reference_method_url' is a string containing the name of
    #' the column in 'table' that contains the reference method urls. Default
    #' column name is 'method_url'.
    #'
    #' In table_names, 'reference_documentation_url' is a string containing the
    #' name of the column in 'table' that contains the reference urls. Default
    #' column name is 'documentation_url'.
    #'
    #' In table_names, 'reference_parameter' is a string containing the name of
    #' the column in 'table' that contains the reference parameters. Default
    #' column name is 'parameter'.
    #'
    #' In table_names, 'reference_citation' is a string containing the name of
    #' the column in 'table' that contains the reference citations. Default
    #' column name is 'citation'.
    #'
    #' @param reference a list created by the function [new_reference].
    #'   Optional.
    #'
    #' @param verbose boolean indicating whether or not you want progress
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
    #' data$add(table = fasta_data, type = "sequence")
    #' data$add(
    #'   table = contigs_report, type = "report",
    #'   report_type = "contigs_report", list(sequence_name = "Name")
    #' )
    #'
    #' # To add metadata related to your study
    #'
    #' metadata <- readRDS(strollur_example("miseq_metadata.rds"))
    #'
    #' data$add(table = metadata, type = "report", report_type = "metadata")
    #'
    #' @return Updated `strollur` object - invisible(self)
    add = function(table,
                   type = "sequence",
                   report_type = NULL,
                   table_names = list(
                     sequence_name = "sequence_name",
                     sequence = "sequence",
                     comment = "comment",
                     reference_vendor = "vendor",
                     reference_name = "name",
                     reference_version = "version",
                     reference_usage = "usage",
                     reference_note = "note",
                     reference_method_url = "method_url",
                     reference_documentation_url = "documentation_url",
                     reference_parameter = "parameter",
                     reference_citation = "citation"
                   ),
                   reference = NULL,
                   verbose = TRUE) {
      default_tn <- list(
        sequence_name = "sequence_name",
        sequence = "sequence",
        comment = "comment",
        reference_vendor = "vendor",
        reference_name = "name",
        reference_version = "version",
        reference_usage = "usage",
        reference_note = "note",
        reference_method_url = "method_url",
        reference_documentation_url = "documentation_url",
        reference_parameter = "parameter",
        reference_citation = "citation"
      )

      table_names <- modifyList(default_tn, table_names)

      # allow for type and report_type to be entered without ""
      type <- as.character(substitute(type))
      if (!is.null(report_type)) {
        report_type <- as.character(substitute(report_type))
      }

      num_added <- 0
      if (type == "sequence") {
        xdev_add_sequences(
          data = self, table = table,
          sequence_name = table_names[["sequence_name"]],
          sequence = table_names[["sequence"]],
          comment = table_names[["comment"]],
          reference = reference,
          verbose = verbose
        )
      } else if (type == "report") {
        if (!is.null(report_type)) {
            # check for sequence name column in table
            if (table_names[["sequence_name"]] %in% base::names(table)) {
                xdev_add_report(
                    data = self, table = table,
                    type = report_type,
                    sequence_name = table_names[["sequence_name"]],
                    verbose = verbose
                )
            } else {
                xdev_add_report(
                    data = self, table = table,
                    type = report_type,
                    verbose = verbose
                )
            }
        } else {
          cli::cli_abort("'report_type' is required when adding a report.")
        }
      } else if (type == "metadata") {
        xdev_add_report(
          data = self, table = table,
          type = type,
          verbose = verbose
        )
      } else if (type == "resource_reference") {
        xdev_add_references(
          data = self, table = table,
          name = table_names[["reference_name"]],
          vendor = table_names[["reference_vendor"]],
          version = table_names[["reference_version"]],
          usage = table_names[["reference_usage"]],
          note = table_names[["reference_note"]],
          method_url = table_names[["reference_method_url"]],
          documentation_url = table_names[["reference_documentation_url"]],
          parameter = table_names[["reference_parameter"]],
          citation = table_names[["reference_citation"]],
          verbose = verbose
        )
      } else {
        message <- paste0(
          type, " is not a valid 'type' for the strollur$add()",
          " function."
        )
        cli::cli_abort(message)
      }

      invisible(self)
    },

    #' @description
    #' Add phylo tree relating the samples in your dataset
    #'
    #' @param tree a phylo tree object created by ape::read.tree.
    #' @examples
    #'
    #'  data <- new_dataset("my_dataset")
    #'
    #'  df <- read_mothur_shared(strollur_example("final.opti_mcc.shared"))
    #'  assign(data = data, table = df, type = "bin", bin_type = "otu")
    #'
    #'  tree <- ape::read.tree(strollur_example(
    #'  "final.opti_mcc.jclass.ave.tre"))
    #'
    #'  data$add_sample_tree(tree)
    #' @return Updated `strollur` object
    add_sample_tree = function(tree) {
      if (!inherits(tree, "phylo")) {
        .abort_incorrect_type("phylo", tree)
      }

      # if no samples, add sequences in tree to dataset
      if (count(self, type = "sample") == 0) {
        message <- paste0("[Warning]: Your dataset does not contain sample ",
          "data, ignoring sample tree.",
          collapse = ""
        )
        cli::cli_alert(message)
      } else {
        # make sure the tree includes all "good" samples
        if (identical(
          sort(tree$tip.label),
          sort(names(self, type = "sample"))
        )) {
          # save tree
          self$sample_tree <- tree
        } else {
          # samples in dataset and not in tree
          missing_samples <- setdiff(
            names(self, type = "sample"),
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
              names(self, type = "sample")
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
    #'  data <- new_dataset("my_dataset")
    #'  tree <- ape::read.tree(strollur_example("final.phylip.tre.gz"))
    #'  data$add_sequence_tree(tree)
    #'
    #' @return Updated `strollur` object
    add_sequence_tree = function(tree) {
      if (!inherits(tree, "phylo")) {
        .abort_incorrect_type("phylo", tree)
      }

      # if no seqs yet, add sequences in tree to dataset
      if (count(self, type = "sequence") == 0) {
        xdev_add_sequences(self, data.frame(sequence_name = tree$tip.label))

        # save tree
        self$sequence_tree <- tree
      } else {
        # make sure the tree includes all "good" sequences
        if (identical(
          sort(tree$tip.label),
          sort(names(self, type = "sequence"))
        )) {
          # save tree
          self$sequence_tree <- tree
        } else {
          # seqs in dataset and not in tree
          missing_seqs <- setdiff(
            names(self, type = "sequence"),
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
              names(self, type = "sequence")
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
    #' @param table a data.frame containing the data you wish to assign
    #'
    #' @param type a string containing the type of data. Options include:
    #' 'sequence_abundance', 'sequence_taxonomy', 'bin',
    #'  'bin_representative', 'bin_taxonomy' and 'treatment'.
    #'  Default = "bin".
    #'
    #' @param bin_type string containing the bin type you would like the number
    #'   of bins for. Default = "otu".
    #'
    #' @param table_names named list used to indicate the names of the columns
    #'   in the table. By default:
    #'
    #'   table_names <- list(sequence_name = "sequence_name", abundance =
    #'   "abundance", sample = "sample", treatment = "treatment", taxonomy =
    #'   "taxonomy", bin_name = "bin_name")
    #'
    #'   In table_names, 'sequence_name' is a string containing the name of the
    #'   column in 'table' that contains the sequence names. Default column name
    #'   is 'sequence_name'.
    #'
    #'   In table_names, 'abundance' is a string containing the name of the
    #'   column in 'table' that contains the abundances. Default column name is
    #'   'abundance'.
    #'
    #'   In table_names, 'sample' is a string containing the name of the column
    #'   in 'table' that contains the samples. Default column name is 'sample'.
    #'
    #'   In table_names, 'treatment' is a string containing the name of the
    #'   column in 'table' that contains the treatment names. Default column
    #'   name is 'treatment'.
    #'
    #'   In table_names, 'taxonomy' is a string containing the name of the
    #'   column in 'table' that contains the classifications. Default column
    #'   name is 'taxonomy'.
    #'
    #'   In table_names, 'bin_name' is a string containing the name of the
    #'   column in 'table' that contains the bin names. Default column name is
    #'   'bin_name'.
    #'
    #' @param reference a list created by the function [new_reference].
    #'   Optional.
    #' @param verbose boolean indicating whether or not you want progress
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
    #' data$assign(table = bin_reps, type = "bin_representative")
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
    #' data$assign(table = sample_assignments, type = "treatment")
    #'
    #' @return Updated `strollur` object
    assign = function(table,
                      type = "bin",
                      bin_type = "otu",
                      table_names = list(
                        sequence_name = "sequence_name",
                        abundance = "abundance",
                        sample = "sample",
                        treatment = "treatment",
                        taxonomy = "taxonomy",
                        bin_name = "bin_name"
                      ),
                      reference = NULL,
                      verbose = TRUE) {
      default_tn <- list(
        sequence_name = "sequence_name",
        abundance = "abundance",
        sample = "sample",
        treatment = "treatment",
        taxonomy = "taxonomy",
        bin_name = "bin_name"
      )

      table_names <- modifyList(default_tn, table_names)

      # allow for type and bin_type to be entered without ""
      type <- as.character(substitute(type))
      bin_type <- as.character(substitute(bin_type))

      num <- 0
      if (type == "bin") {
        num <- xdev_assign_bins(
          data = self, table = table,
          bin_type = bin_type,
          reference = reference,
          bin_name = table_names[["bin_name"]],
          abundance = table_names[["abundance"]],
          sample = table_names[["sample"]],
          sequence_name = table_names[["sequence_name"]],
          verbose = verbose
        )
      } else if (type == "bin_taxonomy") {
        num <- xdev_assign_bin_taxonomy(
          data = self, table = table,
          bin_type = bin_type,
          reference = reference,
          bin_name = table_names[["bin_name"]],
          taxonomy = table_names[["taxonomy"]],
          verbose = verbose
        )
      } else if (type == "bin_representative") {
        num <- xdev_assign_bin_representative_sequences(
          data = self, table = table,
          bin_type = bin_type,
          reference = reference,
          bin_name = table_names[["bin_name"]],
          sequence_name = table_names[["sequence_name"]],
          verbose = verbose
        )
      } else if (type == "sequence_taxonomy") {
        num <- xdev_assign_sequence_taxonomy(
          data = self, table = table,
          reference = reference,
          sequence_name = table_names[["sequence_name"]],
          taxonomy = table_names[["taxonomy"]],
          verbose = verbose
        )
      } else if (type == "treatment") {
        num <- xdev_assign_treatments(
          data = self, table = table,
          sample = table_names[["sample"]],
          treatment = table_names[["treatment"]],
          verbose = verbose
        )
      } else if (type == "sequence_abundance") {
        num <- xdev_assign_sequence_abundance(
          data = self, table = table,
          sequence_name = table_names[["sequence_name"]],
          abundance = table_names[["abundance"]],
          sample = table_names[["sample"]],
          treatment = table_names[["treatment"]],
          verbose = verbose
        )
      } else {
        message <- paste0(
          type, " is not a valid 'type' for the assign()",
          " function."
        )
        cli::cli_abort(message)
      }
      invisible(self)
    },

    #' @description
    #' Clear data from datasest
    #' @examples
    #' miseq <- load_dataset(strollur_example("miseq_sop.rds"))
    #' miseq
    #' miseq$clear()
    #' miseq
    #' @return Updated `strollur` object
    clear = function() {
      clear(self)

      invisible(self)
    },

    #' @description
    #' Find the number of sequences, samples, treatments or bins of a given type
    #'
    #' @param type string containing the type of data you want the number of.
    #' Options include: "sequence", "sample", "treatment", "bin", and
    #'  "resource_reference". Default = "sequence".
    #'
    #' @param bin_type string containing the bin type you would like the number
    #'   of bins for. Default = "otu".
    #'
    #' @param samples vector of strings. samples is only used when 'type' =
    #'   "sequence" or 'type' = "bin" . samples should contain the names of
    #'   the samples you want the count for. Default = NULL.
    #'
    #' @param distinct Boolean. distinct is used when 'type' = "sequence" or
    #'   'type' = "bin". When 'type' = "sequence" and distinct is TRUE the
    #'   number of unique sequences is returned. When 'type' = "sequence" and
    #'   distinct is FALSE the total number of sequences is returned. This can
    #'   also be combined with samples to find the number of unique sequences
    #'   found ONLY in a given set of samples, or to find the number of unique
    #'   sequences in given set of samples that may also be present in other
    #'   samples. When 'type' = "bin", you can set distinct = TRUE to return
    #'   the number of bins that ONLY contain sequences from the given samples.
    #'   When distinct is FALSE the count returned contains bins with sequences
    #'   from a given samples, but those bins may also contain other samples.
    #'   Default = FALSE.
    #'
    #' @examples
    #'
    #' miseq <- load_dataset(strollur_example("miseq_sop.rds"))
    #'
    #' # To get the total number of sequences
    #' miseq$count(type = "sequence")
    #'
    #' # To get number of unique sequences
    #' miseq$count(type = "sequence", distinct = TRUE)
    #'
    #' # To get number of unique sequences from samples 'F3D0' and 'F3D1'
    #' # Note these sequences will be present in both samples but may be
    #' # be present in other samples as well
    #' miseq$count(type = "sequence", samples = c("F3D0", "F3D1"))
    #'
    #' # To get number of unique sequences exclusive to samples 'F3D0' and
    #' # 'F3D1'. Note sequences are present in both samples and NOT present in
    #' # any other samples.
    #'
    #' miseq$count(type = "sequence",
    #'             samples = c("F3D0", "F3D1"), distinct = TRUE )
    #'
    #' # To get the number of samples in the dataset
    #' miseq$count(type = "sample")
    #'
    #' # To get the number of treatments in the dataset
    #' miseq$count(type = "treatment")
    #'
    #' # To get the number of "otu" bins in the dataset
    #' miseq$count(type = "bin", bin_type = "otu")
    #'
    #' # To get the number of "asv" bins in the dataset
    #' miseq$count(type = "bin", bin_type = "asv")
    #'
    #' # To get the number of "phylotype" bins in the dataset
    #' miseq$count(type = "bin", bin_type = "phylotype")
    #'
    #' # To get number of "otu" bins from samples 'F3D0' and 'F3D1'
    #' # Note these bins will have sequences from both samples but there may be
    #' # other samples present as well
    #' miseq$count(
    #'   type = "bin", bin_type = "otu", samples = c("F3D0", "F3D1")
    #' )
    #'
    #' # To get number of "otu" bins unique to samples 'F3D0' and 'F3D1'
    #' # Note these bins will have sequences from both samples and NO other
    #' # samples will be present in the bins.
    #'
    #' miseq$count(
    #'   type = "bin", bin_type = "otu",
    #'   samples = c("F3D0", "F3D1"), distinct = TRUE
    #' )
    #'
    #' @return double
    count = function(type = "sequence",
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
    #'  data <- new_dataset("my_dataset")
    #'
    #'  # assign abundance 'otu' bins
    #'  data$assign(table = df, type = "bin", bin_type = "otu")
    #'
    #'  data$add_sample_tree(tree)
    #'  data$get_sample_tree()
    #'
    #' @return ape::tree
    get_sample_tree = function() {
      if (!is.null(self$sample_tree)) {
        # prune tree if needed
        # samples in tree and not in dataset
        extra_samples <- setdiff(
          self$sample_tree$tip.label,
          names(self, type = "sample")
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
    #'  data <- new_dataset("my_dataset")
    #'  tree <- ape::read.tree(strollur_example("final.phylip.tre.gz"))
    #'  data$add_sequence_tree(tree)
    #'  data$get_sequence_tree()
    #'
    #' @return ape::tree
    get_sequence_tree = function() {
      if (!is.null(self$sequence_tree)) {
        # prune tree if needed
        # seqs in tree and not in dataset
        extra_seqs <- setdiff(
          self$sequence_tree$tip.label,
          names(self, type = "sequence")
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
    #' Get the version of the
    #' \href{https://mothur.org/strollur/reference/strollur.html}{strollur}
    #' object.
    #'
    #' @examples
    #'
    #' data <- new_dataset("test")
    #'
    #' data$get_version()
    #'
    #' @returns a logical
    #' @export
    get_version = function() {
      private$version
    },

    #' @description
    #' Determine if two
    #' \href{https://mothur.org/strollur/reference/strollur.html}{strollur}
    #'  objects are equal.
    #'
    #' @param data, a
    #'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur}
    #'   object
    #' @examples
    #'
    #' miseq <- load_dataset(strollur_example("miseq_sop.rds"))
    #'
    #' data <- copy_dataset(miseq)
    #'
    #' miseq$is_equal(data)
    #'
    #' @returns a logical
    #' @export
    is_equal = function(data) {
      if (!inherits(data, "strollur")) {
        stop("data must be a strollur object.")
      }

      # compare private members - version
      if (!identical(self$.private$version, data$.private$version)) {
        cli_alert("The strollur objects versions are not equivalent.")
        return(FALSE)
      }

      # Compare public members - data, raw, sequence_tree and sample_tree
      # compare raw
      if (!identical(self$raw, data$raw)) {
        message <- paste(
          "The strollur objects public field 'raw' are not",
          " equivalent."
        )
        cli_alert(message)
        return(FALSE)
      }

      # compare ape sequence tree, sample_tree
      # We use ape::all.equal.phylo because tree objects
      # can be topologically identical but have different memory layouts.
      seq_tree_self <- self$get_sequence_tree()
      seq_tree_data <- data$get_sequence_tree()
      if (is.null(seq_tree_self) && is.null(seq_tree_data)) {
        # null trees
      } else if (is.null(seq_tree_self) || is.null(seq_tree_data)) {
        message <- paste(
          "The strollur objects public field ",
          "'sequence_tree' are not equivalent."
        )
        cli_alert(message)
        return(FALSE)
      } else {
        # all.equal.phylo returns TRUE or a character vector of differences
        if (!isTRUE(ape::all.equal.phylo(
          seq_tree_self,
          seq_tree_data
        ))) {
          message <- paste(
            "The strollur objects public field ",
            "'sequence_tree' are not equivalent."
          )
          cli_alert(message)
          return(FALSE)
        }
      }

      sample_tree_self <- self$get_sample_tree()
      sample_tree_data <- data$get_sample_tree()
      if (is.null(sample_tree_self) && is.null(sample_tree_data)) {
        # null trees
      } else if (is.null(sample_tree_self) || is.null(sample_tree_data)) {
        message <- paste(
          "The strollur objects public field ",
          "'sample_tree' are not equivalent."
        )
        cli_alert(message)
        return(FALSE)
      } else {
        # all.equal.phylo returns TRUE or a character vector of differences
        if (!isTRUE(ape::all.equal.phylo(
          sample_tree_self,
          sample_tree_data
        ))) {
          message <- paste(
            "The strollur objects public field ",
            "'sample_tree' are not equivalent."
          )
          cli_alert(message)
          return(FALSE)
        }
      }

      # compare data
      if (!xint_is_equal(self, data)) {
        return(FALSE)
      }

      TRUE
    },


    #' @description
    #' Get the names of a given type of data
    #'
    #' @param type string containing the type of data you would like. Options
    #'   include: "dataset", "sequence", "bin", "sample", "treatment",
    #'   "report". Default = "sequence".
    #'
    #' @param bin_type string containing the bin type you would like the names
    #'   for. Default = "otu".
    #'
    #' @param samples vector of strings. samples is only used when 'type' =
    #'   "sequence" or 'type' = "bin" . samples should contain the names of
    #'   the samples you want names for. Default = NULL.
    #'
    #' @param distinct Boolean. distinct is used when 'type' = "sequence" or
    #'   'type' = "bin" and the samples parameter is used. The distinct
    #'   parameter allows you to get the names that present given set of
    #'   samples. When distinct is TRUE, the names function will return the
    #'   names that ONLY contain data from the given samples. When distinct is
    #'   FALSE the data returned contains data from a given samples, but may
    #'   ALSO contain data from other samples. Default = FALSE.
    #'
    #' @examples
    #'
    #' miseq <- load_dataset(strollur_example("miseq_sop.rds"))
    #'
    #' # To get the name of the dataset
    #' miseq$names(type = "dataset")
    #'
    #' # To get the names of the sequences
    #' miseq$names(type = "sequence")
    #'
    #' # To get the names of the sequences present sample 'F3D0'
    #' miseq$names(type = "sequence", samples = c("F3D0"))
    #'
    #' #' # To get the names of the sequences unique to sample 'F3D0'
    #' miseq$names(type = "sequence", samples = c("F3D0"), distinct = TRUE)
    #'
    #' # To get the names of the samples
    #' miseq$names(type = "sample")
    #'
    #' # To get the names of the treatments
    #' miseq$names(type = "treatment")
    #'
    #' # To get the names of the bins
    #' miseq$names(type = "bin")
    #'
    #' # To get the names of the bins that are unique to 'F3D0'
    #' miseq$names(type = "bin", samples = c("F3D0"), distinct = TRUE)
    #'
    #' # To get the names of the bins that include sequences from 'F3D0'
    #' miseq$names(type = "bin", samples = c("F3D0"), distinct = FALSE)
    #'
    #' # To get the names of the reports
    #' miseq$names(type = "report")
    #'
    #' @return vector of strings, containing the names requested
    names = function(type = "sequence",
                     bin_type = "otu",
                     samples = NULL,
                     distinct = FALSE) {
      xdev_names(self, type, bin_type, samples, distinct)
    },

    #' @description
    #' Get a data.frame containing the given report
    #'
    #' @param type string containing the type of report you would like. Options
    #' include: "fasta", "sequence", "sequence_bin_assignment",
    #' "sequence_taxonomy", "bin_taxonomy", "bin_representative",
    #'  "sample_assignment", "resource_reference", "sequence_scrap",
    #' "bin_scrap". If you have added custom reports for alignment,
    #' contigs_assembly or chimeras, you can get those as well.
    #'  Default = "sequence".
    #'
    #' @param bin_type string containing the bin type you would like a
    #'   bin_taxonomy report for. Default = "otu".
    #'
    #' @examples
    #'
    #' miseq <- load_dataset(strollur_example("miseq_sop.rds"))
    #'
    #' # To get the FASTA data
    #'
    #' miseq$report(type = "fasta") |> head(n = 5)
    #'
    #' # To get a report about the FASTA data
    #'
    #' miseq$report(type = "sequence") |> head(n = 5)
    #'
    #' # To get the sequence bin assignments
    #'
    #' miseq$report(type = "sequence_bin_assignment", bin_type = "otu") |>
    #' head(n = 5)
    #'
    #' # To get the sample treatment assignments
    #'
    #' miseq$report(type = "sample_assignment") |> head(n = 5)
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
    #' miseq$report(type = "bin_representative", bin_type = "otu") |>
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
    #' references <- miseq$report(type = "resource_reference")
    #'
    #' # To get our custom report containing the contigs assembly data:
    #'
    #' miseq$report(type = "contigs_report") |> head(n = 5)
    #'
    #' @return data.frame
    report = function(type = "sequence", bin_type = "otu") {
      xdev_report(self, type, bin_type)
    },

    #' @description
    #' Summarize the sequences data, custom reports, and scrapped data
    #'
    #' @param type string containing the type of data you want the number of.
    #' Options include: "sequence", "report" and "scrap".
    #'  Default = "sequence".
    #'
    #' @param report_type string containing the report type you would
    #'  summarized. For example, the miseq_sop_example includes contigs assembly
    #'  data and can be accessed with report_type = "contigs_report".
    #'  Default = NULL.
    #'
    #' @param verbose boolean indicating whether or not you want progress
    #'  messages. Default = TRUE.
    #'
    #' @examples
    #'
    #' miseq <- load_dataset(strollur_example("miseq_sop.rds"))
    #'
    #' # To get the summary of your FASTA data
    #' miseq$summary(type = "sequence")
    #'
    #' # summarize contigs_report
    #' miseq$summary(type = "report", report_type = "contigs_report")
    #'
    #' # remove sample 'F3D0' to produce a scrap report
    #' xdev_remove_samples(data = miseq, samples = c("F3D0"))
    #'
    #' # summarize scrapped data -
    #' # sequences and bins scrapped by removing the sample "F3D0"
    #' miseq$summary(type = "scrap")
    #'
    #' @return data.frame
    summary = function(type = "sequence",
                       report_type = NULL, verbose = TRUE) {
      dataset_summary <- do.call(summary, list(
        data = self,
        type = type,
        report_type = report_type
      ))
      if (verbose) {
        print(dataset_summary)
      }
      dataset_summary
    }
  ),
  private = list(
    version = "0.1.0",
    finalize = function() {},


    #' Get summary of the sequence reports
    get_summary = function() {
      results <- list()

      # if you have summary results to print
      if (!all(xdev_get_sequences(self) == "")) {
        # if you have summary results to print
        results[["sequence_summary"]] <- summary(
          data = self,
          type = "sequence", verbose = FALSE
        )
      }

      df <- summary(
        data = self,
        type = "scrap", verbose = FALSE
      )
      if (nrow(df) != 0) {
        results[["scrap_summary"]] <- df
      }

      results
    }
  )
)
