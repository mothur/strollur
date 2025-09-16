#' @title sequence_data
#' @description 'sequence_data' is an R6 class that stores nucleotide sequences,
#' abundance, sample and treatment assignments, taxonomic classifications,
#' asv / otu clusters and various reports. It is designed to facilitate data
#' analysis across multiple R packages.
#'
#'
#' @author Sarah Westcott, \email{swestcot@@umich.edu}
#'
#' @importFrom R6 R6Class
#' @importFrom methods new
#' @importFrom parallelly, availableCores
#' @importFrom waldo compare
#' @import cli
#' @export
sequence_data <- R6Class("sequence_data",
  public = list(

    #' @field data Rcpp::XPtr<Dataset> pointer to 'Dataset' c++ class. This
    #' allows package developers an easy access point to the underlying C++ code
    #' with additional functionality.
    data = NULL,

    #' @field raw Rcpp::RawVector containing the serialized data of the
    #' 'Dataset' c++ class. This allows the load and save functions to work with
    #' the class.
    raw = NULL,

    #' @description
    #' Create a new sequence dataset
    #' @param name String, name of dataset (optional)
    #' @param processors Integer, number of cores to use.
    #'  Default = all available
    #' @param dataset a `sequence_data` object.
    #' @examples
    #'
    #' # to create an empty dataset, run the following:
    #'
    #' dataset <- sequence_data$new("soil")
    #'
    #' @return A new `sequence_data` object.
    initialize = function(name = "",
                          processors = parallelly::availableCores(),
                          dataset = NULL) {
      if (is.null(dataset)) {
        self$data <- new_dataset(name, processors)
        private$references <- data.frame()
        private$metadata <- data.frame()
        private$alignment_data <- data.frame()
        private$contigs_data <- data.frame()
        private$sequence_tree <- NULL
        private$processors <- processors
        private$sample_tree <- NULL
      } else {
        # copy of dataset backend
        self$data <- copy_dataset(dataset$data)
        set_num_processors(self$data, processors)
        # assign new name
        if (name != "") {
          set_dataset_name(self$data, name)
        }

        # copy metadata
        private$metadata <- dataset$get_metadata()
        private$references <- dataset$get_references()
        private$alignment_data <- dataset$get_alignment_report()
        private$contigs_data <- dataset$get_contigs_assembly_report()
        private$processors <- processors
        private$sequence_tree <- dataset$get_sequence_tree()
        private$sample_tree <- dataset$get_sample_tree()
      }

      invisible(self)
    },

    #' @description
    #' Get summary of sequence data
    print = function() {
      if (get_dataset_name(self$data) != "") {
        cat(get_dataset_name(self$data))
        cat(":\n\n")
      }
      self$get_sequence_summary()
      cat("\n\n")
      self$get_sample_summary()
      if (self$get_num_sequences(TRUE) != 0) {
        cat(
          paste("\nNumber of unique seqs:", self$get_num_sequences(TRUE)),
          "\n"
        )
      } else {
        cat("\n")
      }
      cat(
        paste("Total number of seqs:", self$get_num_sequences()),
        "\n"
      )

      if (self$get_num_bins("otu") != 0) {
        cat(
          paste("Total number of otus:", self$get_num_bins("otu")),
          "\n"
        )
      }
      if (self$get_num_bins("asv") != 0) {
        cat(
          paste(
            "Total number of asvs:", self$get_num_bins("asv"),
            "\n"
          )
        )
      }
      if (self$get_num_bins("phylotype") != 0) {
        cat(
          paste(
            "Total number of phylotype bins:",
            self$get_num_bins("phylotype"),
            "\n"
          )
        )
      }
      cat("\n")
    },

    #' @description
    #' Add alignment report for your dataset
    #' @param report a data.frame containing alignment data about your dataset
    #' @param sequence_name_column a string containing the name of the column
    #' containing the sequence names.
    #' @examples
    #'
    #'   dataset <- sequence_data$new("my_dataset")
    #'   align_report <- readr::read_tsv(rdataset_example("alignment_data.tsv"),
    #'    col_names = TRUE, show_col_types = FALSE)
    #'   dataset$add_alignment_report(align_report, "QueryName")
    #'
    add_alignment_report = function(report, sequence_name_column) {
      if (!is.data.frame(report)) {
        abort_incorrect_type("data.frame", report)
      }

      if (!(sequence_name_column %in% names(report))) {
        message <- paste("[ERROR]: The alignment report must include a ",
          "column containing sequence names. ",
          sequence_name_column,
          " is not a named column in your alignment_report.",
          collapse = ""
        )
        cli::cli_abort(message)
      } else {
        # no sequence data yet, add the sequences
        if (self$get_num_sequences() == 0) {
          self$add_sequences(sequence_names = report[[sequence_name_column]])
          private$alignment_data <- report
          # save name column
          attr(
            private$alignment_data,
            "sequence_name_column"
          ) <- sequence_name_column
        } else {
          # seqs in dataset and not in report
          missing_seqs <- setdiff(
            self$get_sequence_names(),
            report[[sequence_name_column]]
          )

          # make sure there is a report entry for each sequence in dataset
          if (length(missing_seqs) == 0) {
            # preserve order of dataset
            private$alignment_data <- report[order(
              match(
                report[[sequence_name_column]],
                self$get_sequence_names()
              )
            ), ]
            # save name column
            attr(
              private$alignment_data,
              "sequence_name_column"
            ) <- sequence_name_column
          } else {
            message <- paste("[WARNING]: Your alignment report does no",
              "t contain an entry for every sequence in",
              " your dataset, ignoring alignment report",
              ". Missing alignment report entries for: ",
              paste(missing_seqs, collapse = ", "),
              ".",
              collapse = ""
            )
            cli_alert(message)
          }
        }
      }

      invisible(self)
    },

    #' @description
    #' Add contigs assembly report for your dataset
    #' @param report a data.frame containing contigs assembly data about your
    #' dataset
    #' @param sequence_name_column a string containing the name of the column
    #' containing the sequence names.
    #' @examples
    #'
    #'   dataset <- sequence_data$new("my_dataset")
    #'   contigs_report <- readr::read_tsv(rdataset_example("contigs_data.tsv"),
    #'    col_names = TRUE, show_col_types = FALSE)
    #'   dataset$add_contigs_assembly_report(contigs_report, "Name")
    #'
    add_contigs_assembly_report = function(report, sequence_name_column) {
      if (!is.data.frame(report)) {
        abort_incorrect_type("data.frame", report)
      }

      if (!(sequence_name_column %in% names(report))) {
        message <- paste("[ERROR]: The contigs assembly report must ",
          "include a column containing sequence names. ",
          sequence_name_column,
          " is not a named column in your contigs_report.",
          collapse = ""
        )
        cli::cli_abort(message)
      } else {
        # no sequence data yet, add the sequences
        if (self$get_num_sequences() == 0) {
          self$add_sequences(sequence_names = report[[sequence_name_column]])
          private$contigs_data <- report
          # save name column
          attr(
            private$contigs_data,
            "sequence_name_column"
          ) <- sequence_name_column
        } else {
          # seqs in dataset and not in report
          missing_seqs <- setdiff(
            self$get_sequence_names(),
            report[[sequence_name_column]]
          )

          # make sure there is a report entry for each sequence in dataset
          if (length(missing_seqs) == 0) {
            # preserve order of dataset
            private$contigs_data <- report[order(
              match(
                report[[sequence_name_column]],
                self$get_sequence_names()
              )
            ), ]
            # save name column
            attr(
              private$contigs_data,
              "sequence_name_column"
            ) <- sequence_name_column
          } else {
            message <- paste("[WARNING]: Your contigs assembly report ",
              "does not contain an entry for every ",
              "sequence in your dataset, ignoring. ",
              "Missing contigs",
              " assembly report entries for: ",
              paste(missing_seqs, collapse = ", "),
              ".",
              collapse = ""
            )
            cli_alert(message)
          }
        }
      }

      invisible(self)
    },

    #' @description
    #' Add metadata about your dataset
    #' @param metadata a data.frame containing metadata about your dataset
    #' @examples
    #'
    #'   dataset <- sequence_data$new("my_dataset")
    #'   metadata <- readr::read_tsv(rdataset_example("sample-metadata.tsv"),
    #'    col_names = TRUE, show_col_types = FALSE)
    #'   dataset$add_metadata(metadata)
    #'
    add_metadata = function(metadata) {
      if (!is.data.frame(metadata)) {
        abort_incorrect_type("data.frame", metadata)
      }

      private$metadata <- metadata

      invisible(self)
    },

    #' @description
    #' Add information about the references used to analyze your dataset
    #' @param reference a data.frame containing multiple references. The only
    #' required column is reference_name. Optional columns can include: version,
    #' usage, note and url. Default = NULL. Either the reference or
    #' reference_name parameter is required.
    #' @param reference_name a string containing the name of the reference. For
    #' example: 'silva.bacteria.fasta' Default = NULL. Either the reference
    #' or reference_name parameter is required.
    #' @param version a string containing the version of the reference. For
    #' eaxmple: '138.2' Default = NULL.
    #' @param usage a string containing the usage of the reference in
    #' your analysis. For example: 'alignment using mothur2' Default = NULL.
    #' @param note a string containing the any additional notes about the
    #' reference. For example: 'custom reference using silva version 1.38.2 with
    #' additional sequences tailored for this analysis'. Default = NULL.
    #' @param url a string containing a web address where the reference may be
    #' downloaded. For example: 'https://mothur.org/wiki/silva_reference_files/'
    #' . Default = NULL.
    #' @examples
    #'
    #' # To add a single reference:
    #' dataset <- sequence_data$new("my_dataset")
    #' dataset$add_references(reference_name = "silva.v4.fasta",
    #' version = "1.38.1", usage = "alignment", note = "custom reference
    #'  created by trimming silva.bacteria.fasta to the V4 region",
    #'  url = "https://mothur.org/wiki/silva_reference_files/")
    #'
    #' # To add multiple references:
    #' reference <- readr::read_csv(rdataset_example("references.csv"),
    #'                              col_names = TRUE, show_col_types = FALSE)
    #' dataset$add_references(reference = reference)
    #'
    add_references = function(reference = NULL, reference_name = NULL,
                              version = NULL, usage = NULL, note = NULL,
                              url = NULL) {
      # make sure at least one reference is given
      if (is.null(reference) && is.null(reference_name)) {
        parameters <- list("reference", "reference_name")
        abort_provide_at_least_one(parameters)
      }

      # add multiple references
      if (!is.null(reference)) {
        if (!is.data.frame(reference)) {
          abort_incorrect_type("data.frame", reference)
        }

        if (!("reference_name" %in% names(reference))) {
          cli::cli_abort("[ERROR]: The reference data.frame must include
                               a 'reference_name' column.")
        } else {
          rows <- nrow(private$references)
          if (rows == 0) {
            private$references <- reference
          } else {
            # find missing columns
            missing_columns <- setdiff(
              names(reference),
              names(private$references)
            )

            # add missing columns to existing references
            for (col in missing_columns) {
              private$references[rows, col] <- NA
            }

            private$references <- rbind(private$references, reference)
          }
        }
      }

      # add single reference if provided
      if (!is.null(reference_name)) {
        rows <- nrow(private$references) + 1

        # new blank row with all columns in the reference
        private$references[rows, ] <- NA

        if (is.null(version)) {
          version <- NA
        }

        if (is.null(usage)) {
          usage <- NA
        }

        if (is.null(note)) {
          note <- NA
        }

        if (is.null(url)) {
          url <- NA
        }

        private$references[rows, "reference_name"] <- reference_name
        private$references[rows, "version"] <- version
        private$references[rows, "usage"] <- usage
        private$references[rows, "note"] <- note
        private$references[rows, "url"] <- url
      }

      invisible(self)
    },

    #' @description
    #' Add phylo tree relating the samples in your dataset
    #' @param tree a phylo tree object created by ape::read.tree.
    #' @examples
    #'
    #'  dataset <- sequence_data$new("my_dataset")
    #'  df <- read_mothur_shared(rdataset_example("final.opti_mcc.shared"))
    #'  dataset$assign_bins(df$bin_id, df$abundance, df$sample)
    #'  tree <- ape::read.tree(rdataset_example(
    #'  "final.opti_mcc.jclass.ave.tre"))
    #'  dataset$add_sample_tree(tree)
    #'
    add_sample_tree = function(tree) {
      if (!inherits(tree, "phylo")) {
        abort_incorrect_type("phylo", tree)
      }

      # if no samples, add sequences in tree to dataset
      if (self$get_num_samples() == 0) {
        message <- paste0("[Warning]: Your dataset does not contain sample ",
          "data, ignoring sample tree.",
          collapse = ""
        )
        cli::cli_alert(message)
      } else {
        # make sure the tree includes all "good" samples
        if (identical(
          sort(tree$tip.label),
          sort(self$get_samples())
        )) {
          # save tree
          private$sample_tree <- tree
        } else {
          # samples in dataset and not in tree
          missing_samples <- setdiff(
            self$get_samples(),
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
              self$get_samples()
            )

            # if tree contains "extra" names, prune the tree
            private$sample_tree <- drop.tip(tree, tip = extra_samples)
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
    #'  dataset <- sequence_data$new("my_dataset")
    #'  tree <- ape::read.tree(rdataset_example("final.phylip.tre"))
    #'  dataset$add_sequence_tree(tree)
    #'
    add_sequence_tree = function(tree) {
      if (!inherits(tree, "phylo")) {
        abort_incorrect_type("phylo", tree)
      }

      # if no seqs yet, add sequences in tree to dataset
      if (self$get_num_sequences() == 0) {
        self$add_sequences(sequence_names = tree$tip.label)

        # save tree
        private$sequence_tree <- tree
      } else {
        # make sure the tree includes all "good" sequences
        if (identical(
          sort(tree$tip.label),
          sort(self$get_sequence_names())
        )) {
          # save tree
          private$sequence_tree <- tree
        } else {
          # seqs in dataset and not in tree
          missing_seqs <- setdiff(
            self$get_sequence_names(),
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
              self$get_sequence_names()
            )

            # if tree contains "extra" names, prune the tree
            private$sequence_tree <- drop.tip(tree, tip = extra_seqs)
          }
        }
      }

      invisible(self)
    },

    #' @description
    #' Add new sequence data
    #'
    #' @param data a data.frame containing names, sequences(optional) and
    #' comments(optional).
    #' @param sequence_names a vector of sequence names or if using the 'data'
    #' parameter a string containing the name of the column in 'data' that
    #' contains the sequence names. Default column name is 'names'.
    #' @param sequences a vector of sequence nucleotide strings or if using the
    #' 'data' parameter a string containing the name of the column in 'data'
    #' that contains the sequence nucleotide strings. Default column name is
    #'  'sequences'.
    #' @param comments a vector of sequence comments or if using the
    #' 'data' parameter a string containing the name of the column in 'data'
    #' that contains the sequence comments. Default column name is 'comments'.
    #'
    #' @param reference_name a string containing the name of the reference used
    #' in the classification of the sequences. For example:
    #' 'silva.bacteria.fasta' Default = NULL. (optional)
    #' @param reference_version a string containing the version of the reference
    #' used in the preparing of the sequences. For example: '1.38.1'.
    #' Default = NULL. (optional)
    #' @param reference_note a string containing the any additional notes about
    #' the reference. For example: 'used for alignment using mothur2 v1.0'.
    #' Default = NULL. (optional)
    #' @param reference_url a string containing a web address where the
    #' reference may be downloaded. Default = NULL. (optional)
    #' @examples
    #'
    #'   dataset <- sequence_data$new("my_dataset")
    #'   fasta_data <- read_fasta(rdataset_example("final.fasta"))
    #'   dataset$add_sequences(fasta_data)
    #'
    #' # With the additional parameters to add information about the reference
    #' # You can also add references using the 'add_references' function.
    #'
    #' url <- "https://mothur.org/wiki/silva_reference_files/"
    #'
    #' dataset <- sequence_data$new("my_dataset")
    #' dataset$add_sequences(fasta_data,
    #' reference_name = "silva.bacteria.fasta",
    #' reference_note = "alignment by mothur2 v1.0 using default options",
    #' reference_version = "1.38.1", reference_url = url)
    #'
    add_sequences = function(data = NULL,
                             sequence_names = NULL, sequences = NULL,
                             comments = NULL,
                             reference_name = NULL,
                             reference_version = NULL,
                             reference_note = NULL,
                             reference_url = NULL) {
      if (is.null(data) && (is.null(sequence_names))) {
        abort_provide_at_least_one(c("data", "sequence_names"))
      }

      if (!is.null(data)) {
        # required
        sequence_names <- private$fill_required_param(
          sequence_names, data,
          "names"
        )

        # optional
        sequences <- private$fill_optional_param(sequences, data, "sequences")

        # optional
        comments <- private$fill_optional_param(comments, data, "comments")
      } else {
        if (is.null(comments)) {
          comments <- c("")
        }
        if (is.null(sequences)) {
          sequences <- c("")
        }
      }

      add_sequences(self$data, sequence_names, sequences, comments)

      # if a reference is given, save it
      if (!is.null(reference_name)) {
        self$add_references(
          reference_name = reference_name,
          version = reference_version,
          note = reference_note,
          url = reference_url
        )
      }

      invisible(self)
    },

    #' @description
    #' Add bin assignments
    #'
    #' @param data a data.frame containing bin_names, abundances(optional) and
    #' samples(optional), sequence_names(optional). You must provide either
    #' abundances or sequence_names.
    #'
    #' @param bin_names a vector strings containing bin labels or if using the
    #' 'data' parameter a string containing the name of the column in 'data'
    #' that contains the bin names. Default column name is 'bin_names'.
    #' (required)
    #' @param abundances a vector of abundances or if using the 'data' parameter
    #'  a string containing the name of the column in 'data' that contains
    #'  the abundances. Default column name is 'abundances'. You must provide
    #'  either 'abundances' or 'sequence_names'.
    #' @param samples a vector of strings containing sample assignments or if
    #'  using the 'data' parameter a string containing the name of the column
    #'  in 'data' that contains the samples. Default column name is 'samples'.
    #' @param sequence_names a vector of strings containing sequence names or if
    #' using the 'data' parameter a string containing the name of the column in
    #' 'data' that contains the sequence names. Default column name is
    #' 'sequence_names'. You must provide either 'abundances' or
    #' 'sequence_names'.
    #' @param type a string indicating the type of bin assignments.
    #' Default = "otu".
    #' @examples
    #'
    #'   # To assign sequences to bins:
    #'
    #'   dataset <- sequence_data$new("my_dataset")
    #'   seq_ids <- c("seq1", "seq2", "seq4", "seq3", "seq6", "seq5")
    #'   bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
    #'   dataset$assign_bins(bin_names = bin_ids, sequence_names = seq_ids)
    #'
    #'   # bins would look like:
    #'   #            bin1             bin2        bin3
    #'   # (list)     seq1,seq2,seq4   seq3,seq6   seq5
    #'
    #'   # To add abundance only bin assignments:
    #'
    #'   dataset <- sequence_data$new("my_dataset")
    #'   bin_ids <- c("bin1", "bin2", "bin3")
    #'   abundances <- c(110, 525, 80)
    #'   dataset$assign_bins(bin_names = bin_ids, abundances)
    #'
    #'   # bins would look like:
    #'   #            bin1             bin2        bin3
    #'   # (rabund)   110              525         80
    #'
    #'   # To add abundance bin assignments parsed by sample:
    #'
    #'   dataset <- sequence_data$new("my_dataset")
    #'   bin_table <- readr::read_tsv(rdataset_example(
    #'                                "mothur2_bin_assignments_shared.tsv")))
    #'   dataset$assign_bins(bin_table, bin_names = "id",
    #'                       abundances = "abundance", samples = "sample")
    #'
    #'   # To assign sequences to bins with their abundances parsed by sample:
    #'
    #'   dataset <- sequence_data$new("my_dataset")
    #'   bin_ids <- c("bin1", "bin1", "bin1", "bin1", "bin1", "bin1",
    #'                "bin2", "bin2", "bin2",
    #'                "bin3", "bin3")
    #'   seq_ids <- c("seq1", "seq1", "seq1", "seq2", "seq4", "seq4",
    #'                "seq3", "seq3", "seq6",
    #'                "seq5", "seq5")
    #'   samples <- c("sample1", "sample2", "sample5",
    #'                "sample1", "sample3", "sample4",
    #'                "sample2", "sample3", "sample1",
    #'                "sample1", "sample6")
    #'   abundances <- c(10, 100, 1, 500, 25, 80, 20, 5, 60, 15, 50)
    #'   dataset$assign_bins(bin_names = bin_ids, abundances, samples, seq_ids)
    #'
    #'   # bins would look like:
    #'   #            bin1             bin2        bin3
    #'   # (list)     seq1,seq2,seq4   seq3,seq6   seq5
    #'   # (rabund)   716              85          65
    #'
    #'   dataset$get_bin_assignments()
    #'
    assign_bins = function(data = NULL, bin_names = NULL, abundances = NULL,
                           samples = NULL, sequence_names = NULL,
                           type = "otu") {
      if (is.null(data) && (is.null(bin_names))) {
        abort_provide_at_least_one(c("data", "bin_names"))
      }

      if (!is.null(data)) {
        # required
        bin_names <- private$fill_required_param(
          bin_names, data,
          "bin_names"
        )
        # optional
        abundances <- private$fill_optional_param(
          abundances, data,
          "abundances"
        )
        # optional
        samples <- private$fill_optional_param(samples, data, "samples")

        # optional
        sequence_names <- private$fill_optional_param(
          sequence_names, data,
          "sequence_names"
        )

        # neither abundances or sequence_names were provided or found
        lbn <- length(bin_names)
        if ((lbn != length(abundances)) && (lbn != length(sequence_names))) {
          cli::cli_abort("[ERROR]: You must provide either
                                abundances or sequence_names")
        }else if (lbn != length(abundances)) {
            abundances <- 0
        }
      } else {
        if (is.null(abundances) && is.null(sequence_names)) {
          cli::cli_abort("[ERROR]: You must provide either
                                abundances or sequence_names")
        }

        if (is.null(samples)) {
          samples <- ""
        }
        if (is.null(abundances)) {
          abundances <- 0
        }
        if (is.null(sequence_names)) {
          sequence_names <- ""
        }
      }

      assign_bins(
        self$data,
        bin_names, abundances,
        samples, sequence_names, type
      )

      invisible(self)
    },

    #' @description
    #' Assign bin classification.
    #'
    #' Note, if you assign sequence taxonomies and assign bins, 'sequence_data'
    #' will find the concensus taxonomy for each bin for you.
    #' @param bin_ids a vector of bin names
    #' @param taxonomies a vector of bin classifications
    #' @param type a string indicating the type of clusters. Default = "otu".
    #' @param reference_name a string containing the name of the reference used
    #' in the classification of the bins. For example:
    #' 'trainset9_032012.pds.zip' Default = NULL.
    #' @param reference_version a string containing the version of the reference
    #' used in the classification of the bins For example: '9_032012'.
    #' Default = NULL.
    #' @param reference_note a string containing the any additional notes about
    #' the reference. For example: 'custom reference based on RDP'.
    #' Default = NULL.
    #' @param reference_url a string containing a web address where the
    #' reference may be downloaded. Default = NULL.
    #' @examples
    #'
    #' bin_ids <- c("bin1", "bin2", "bin3", "bin4")
    #' abunds <- c(200, 40, 100, 5)
    #' taxonomies <- c("Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;",
    #'               "Bacteria;Proteobacteria;Betaproteobacteria;Neisseriales;",
    #'                "Bacteria;Firmicutes;Bacilli;Lactobacillales;",
    #'            "Bacteria;Proteobacteria;Gammaproteobacteria;Pasteurellales;")
    #'
    #' dataset <- sequence_data$new("my_dataset")
    #' dataset$assign_bins(bin_ids, abunds)
    #' dataset$assign_bin_taxonomy(bin_ids, taxonomies)
    #'
    #' # With the additional parameters to add information about the reference
    #' # You can also add references using the 'add_references' function.
    #'
    #' url <- paste("https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset",
    #'          "9_032012.pds.zip", collapse = "")
    #'
    #' dataset$assign_bin_taxonomy(bin_ids, taxonomies,
    #' reference_name = "trainset9_032012.pds.zip",
    #' reference_note = "classification by mothur2 v1.0 using default options",
    #' reference_version = "9_032012", reference_url = url)
    #'
    assign_bin_taxonomy = function(bin_ids, taxonomies, type = "otu",
                                   reference_name = NULL,
                                   reference_version = NULL,
                                   reference_note = NULL,
                                   reference_url = NULL) {
      assign_bin_taxonomy(self$data, bin_ids, taxonomies, type)

      # if a reference is given, save it
      if (!is.null(reference_name)) {
        self$add_references(
          reference_name = reference_name,
          version = reference_version,
          usage = "bin_classification",
          note = reference_note,
          url = reference_url
        )
      }

      invisible(self)
    },

    #' @description
    #' Add sequence abundance data with optional sample / treatment assignments
    #'
    #' @param data a data.frame containing sequence_names, abundances,
    #' samples(optional), treatments(optional). Either 'data' or
    #'  'sequence_names' and 'abundances' is required.
    #'
    #' @param sequence_names a vector of sequence names or if using the 'data'
    #' parameter a string containing the name of the column in 'data' that
    #' contains the sequence names. Default column name is 'names'. Either
    #' 'data' or 'sequence_names' and 'abundances' is required.
    #' @param abundances a vector of sequence abundances or if using the 'data'
    #' parameter a string containing the name of the column in 'data' that
    #' contains the sequence abundances. Default column name is 'abundances'.
    #' Either 'data' or 'sequence_names' and 'abundances' is required.
    #' @param samples a vector of sample names or if using the 'data'
    #' parameter a string containing the name of the column in 'data' that
    #' contains the sample names. Default column name is 'samples'.
    #' (optional)
    #' @param treatments a vector of treatment assignments or if using the
    #' 'data' parameter a string containing the name of the column in 'data'
    #' that contains the treatments names. Default column name is 'treatments'.
    #' (optional)
    #'
    #' @examples
    #'
    #' dataset <- sequence_data$new("miseq_sop")
    #' sequence_abundance <- readr::read_tsv(rdataset_example(
    #'                                       "mothur2_count_table.tsv"))
    #' dataset$assign_sequence_abundance(sequence_abundance)
    #'
    assign_sequence_abundance = function(data = NULL, sequence_names = NULL,
                                         abundances = NULL, samples = NULL,
                                         treatments = NULL) {
      if (is.null(data) && (is.null(sequence_names))) {
        abort_provide_at_least_one(c("data", "sequence_names"))
      }

      if (!is.null(data)) {
        # required
        sequence_names <- private$fill_required_param(
          sequence_names, data,
          "names"
        )
        # required
        abundances <- private$fill_required_param(
          abundances, data,
          "abundances"
        )
        # optional
        samples <- private$fill_optional_param(samples, data, "samples")

        # optional
        treatments <- private$fill_optional_param(
          treatments, data,
          "treatments"
        )
      } else {
        if (is.null(sequence_names) || (is.null(abundances))) {
          cli::cli_abort("[ERROR]: Unless using the data parameter,
                         'sequence_names' and 'abundances' are required.")
        }

        if (length(sequence_names) != length(abundances)) {
          cli::cli_abort("[ERROR]: sequence_names and abundances must be the
                         same length.")
        }

        if (is.null(samples)) {
          samples <- ""
        }

        if (is.null(treatments)) {
          treatments <- ""
        }
      }

      assign_sequence_abundance(
        self$data, sequence_names, abundances,
        samples, treatments
      )
      invisible(self)
    },

    #' @description
    #' Assign sequence classification
    #' @param names a vector of sequence names
    #' @param taxonomies a vector of sequence classifications
    #' @param reference_name a string containing the name of the reference used
    #' in the classification of the sequences. For example:
    #' 'trainset9_032012.pds.zip' Default = NULL.
    #' @param reference_version a string containing the version of the reference
    #' used in the classification of the sequences. For example: '9_032012'.
    #' Default = NULL.
    #' @param reference_note a string containing the any additional notes about
    #' the reference. For example: 'custom reference based on RDP'.
    #' Default = NULL.
    #' @param reference_url a string containing a web address where the
    #' reference may be downloaded. Default = NULL.
    #' @examples
    #'
    #' names <- c("seq1", "seq2", "seq3", "seq4")
    #' taxonomies <- c("Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;",
    #'               "Bacteria;Proteobacteria;Betaproteobacteria;Neisseriales;",
    #'                "Bacteria;Firmicutes;Bacilli;Lactobacillales;",
    #'            "Bacteria;Proteobacteria;Gammaproteobacteria;Pasteurellales;")
    #'
    #' dataset <- sequence_data$new("my_dataset")
    #' dataset$assign_sequence_taxonomy(names, taxonomies)
    #'
    #' # With the additional parameters to add information about the reference
    #' # You can also add references using the 'add_references' function.
    #'
    #' url <- paste("https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset",
    #'          "9_032012.pds.zip", collapse = "")
    #'
    #' dataset$assign_sequence_taxonomy(names = names, taxonomies = taxonomies,
    #' reference_name = "trainset9_032012.pds.zip",
    #' reference_note = "classification by mothur2 v1.0 using default options",
    #' reference_version = "9_032012", reference_url = url)
    #'
    assign_sequence_taxonomy = function(names, taxonomies,
                                        reference_name = NULL,
                                        reference_version = NULL,
                                        reference_note = NULL,
                                        reference_url = NULL) {
      assign_sequence_taxonomy(self$data, names, taxonomies)

      # if a reference is given, save it
      if (!is.null(reference_name)) {
        self$add_references(
          reference_name = reference_name,
          version = reference_version,
          usage = "sequence_classification",
          note = reference_note,
          url = reference_url
        )
      }

      invisible(self)
    },

    #' @description
    #' Assign samples to treatments
    #' @param samples a vector of sample names
    #' @param treatments a vector of treatment names
    #' @examples
    #'
    #' names <- c("seq1", "seq1", "seq1", "seq2", "seq2",
    #'              "seq2", "seq3", "seq3", "seq4")
    #' samples <- c("sample2", "sample3", "sample4",
    #'             "sample2", "sample3", "sample4",
    #'            "sample2", "sample3",
    #'            "sample4")
    #' abundances <- c(250, 400, 500,
    #'                25, 40, 50,
    #'                25, 25,
    #'                4)
    #' treatments <- c("early", "early", "late")
    #'
    #' dataset <- sequence_data$new("my_dataset")
    #' dataset$assign_sequence_abundance(data = NULL,
    #'                                   names, abundances, samples)
    #' dataset$assign_treatments(unique(samples), treatments)
    #'
    assign_treatments = function(samples, treatments) {
      assign_treatments(self$data, samples, treatments)
      invisible(self)
    },

    #' @description
    #' Remove 'all', 'metadata', 'references', 'alignment_report' or
    #' 'contigs_assembly_report' data from your dataset.
    #' @param type a string indicating the type of data you want to remove.
    #' Options include: "all", "metadata", "references", "alignment_report" or
    #' "contigs_assembly_report". Default = "all".
    clear = function(type = "all") {
      if (type == "all") {
        clear(self$data)
        private$metadata <- data.frame()
        private$references <- data.frame()
        private$alignment_data <- data.frame()
        private$contigs_data <- data.frame()
        private$sequence_tree <- NULL
        private$sample_tree <- NULL
        self$raw <- NULL
      } else if (type == "metadata") {
        private$metadata <- data.frame()
      } else if (type == "references") {
        private$references <- data.frame()
      } else if (type == "alignment_report") {
        private$alignment_data <- data.frame()
      } else if (type == "contigs_assembly_report") {
        private$contigs_data <- data.frame()
      } else {
        cli_alert("{.var {type}} is not a valid type to clear, ignoring.")
      }

      invisible(self)
    },

    #' @description
    #' Get List containing dataset
    #' @return List
    export = function() {
      export(self$data)
    },

    #' @description
    #' Get data.frame containing the alignment report data
    #' @examples
    #'
    #'   dataset <- sequence_data$new("my_dataset")
    #'   align_report <- readr::read_tsv(rdataset_example("alignment_data.tsv"),
    #'    col_names = TRUE, show_col_types = FALSE)
    #'   dataset$add_alignment_report(align_report, "QueryName")
    #'   dataset$get_alignment_report()
    #'
    #' @return data.frame
    get_alignment_report = function() {
      # if there is alignment data
      if (nrow(private$alignment_data) != 0) {
        name_col <- attr(private$alignment_data, "sequence_name_column")

        # select alignment data for the "good" sequences
        df <- private$alignment_data[
          private$alignment_data[[name_col]] %in%
            self$get_sequence_names(),
        ]
        # match order to sequences
        return(as.data.frame(df[order(match(
          df[[name_col]],
          self$get_sequence_names()
        )), ]))
      }
      data.frame()
    },

    #' @description
    #' Get data frame containing sequence bin assignments
    #' @param type a string indicating the type of clusters. Options
    #' include: "otu", "asv", or "phylotype". Default = "otu".
    #' @examples
    #'   dataset <- sequence_data$new("my_dataset")
    #'   bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
    #'   samples <- c("sample1", "sample2", "sample5",
    #'    "sample1", "sample3", "sample1")
    #'   sample_abundances <- c(10, 100, 1, 500, 25, 80)
    #'   dataset$assign_bins(bin_ids, sample_abundances, samples)
    #'
    #'   # (shared) bins would look like:
    #'   # label  sample   bin1   bin2   bin3
    #'   # 0.03   sample1  10     500    80
    #'   # 0.03   sample2  100    0      0
    #'   # 0.03   sample3  0      25     0
    #'   # 0.03   sample5  1      0      0
    #'
    #'   shared <- dataset$get_bin_assignments()
    #'
    #'   # shared$bin_id  shared$sample  shared$abundance
    #'   # bin1           sample1        10
    #'   # bin1           sample2        100
    #'   # bin1           sample5        1
    #'   # bin2           sample1        500
    #'   # bin2           sample3        25
    #'   # bin3           sample1        80
    #'
    #' @return data.frame
    get_bin_assignments = function(type = "otu") {
      get_bin_assignments(self$data, type)
    },

    #' @description
    #' Get report containing the bin taxonomy table -
    #' ids, taxonomy by levels
    #' @param type a string indicating the type of bin clusters.
    #' Default = "otu".
    #' @examples
    #'
    #' bin_ids <- c("bin1", "bin2", "bin3", "bin4")
    #' abunds <- c(200, 40, 100, 5)
    #' taxonomies <- c("Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;",
    #'               "Bacteria;Proteobacteria;Betaproteobacteria;Neisseriales;",
    #'                "Bacteria;Firmicutes;Bacilli;Lactobacillales;",
    #'            "Bacteria;Proteobacteria;Gammaproteobacteria;Pasteurellales;")
    #'
    #' dataset <- sequence_data$new("my_dataset")
    #' dataset$assign_bins(bin_ids, abunds)
    #' dataset$assign_bin_taxonomy(bin_ids, taxonomies)
    #' dataset$get_bin_taxonomy_report()
    #'
    #' @return data.frame
    get_bin_taxonomy_report = function(type = "otu") {
      get_bin_taxonomy_report(self$data, type)
    },

    #' @description
    #' Get data.frame containing contigs report
    #' @examples
    #'
    #'   dataset <- sequence_data$new("my_dataset")
    #'   contigs_report <- readr::read_tsv(rdataset_example("contigs_data.tsv"),
    #'    col_names = TRUE, show_col_types = FALSE)
    #'   dataset$add_contigs_assembly_report(contigs_report, "Name")
    #'   dataset$get_contigs_assembly_report()
    #'
    #' @return data.frame
    get_contigs_assembly_report = function() {
      # if there is contigs data
      if (nrow(private$contigs_data) != 0) {
        name_col <- attr(private$contigs_data, "sequence_name_column")

        # select alignment data for the "good" sequences
        df <- private$contigs_data[
          private$contigs_data[[name_col]] %in%
            self$get_sequence_names(),
        ]
        # match order to sequences
        return(as.data.frame(df[order(match(
          df[[name_col]],
          self$get_sequence_names()
        )), ]))
      }
      data.frame()
    },

    #' @description
    #' Get count table returns data.frame containing:
    #' ids, abundances, sample(optional), treatment(optional)
    #' This table represents mothur's count and design files.
    #' @return data.frame
    get_count_table = function() {
      get_sequence_abundance_table(self$data)
    },

    #' @description
    #' Get dataset name
    #' @return String
    get_dataset_name = function() {
      get_dataset_name(self$data)
    },

    #' @description
    #' Get data frame containing sequence bin assignments
    #' @param type a string indicating the type of clusters. Options
    #' include: "otu", "asv", or "phylotype". Default = "otu".
    #' @examples
    #'   dataset <- sequence_data$new("my_dataset")
    #'   seq_ids <- c("seq1", "seq2", "seq4", "seq3", "seq6", "seq5")
    #'   sequence_abundances <- c(10, 100, 1, 500, 25, 80)
    #'   bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
    #'   dataset$assign_bins(bin_ids, sequence_abundances,
    #'    sequence_names = seq_ids)
    #'
    #'   # (list) bins would look like:
    #'   # bin1             bin2        bin3
    #'   # seq1,seq2,seq4   seq3,seq6   seq5
    #'
    #'   list <- dataset$get_list()
    #'
    #'   #  list$bin_id  list$seq_id
    #'   #  bin1          seq1
    #'   #  bin1          seq2
    #'   #  bin1          seq4
    #'   #  bin2          seq3
    #'   #  bin2          seq6
    #'   #  bin3          seq5
    #'
    #' @return data.frame
    get_list = function(type = "otu") {
      get_list(self$data, type)
    },

    #' @description
    #' Get data.frame containing metadata for the dataset
    #' @examples
    #'   dataset <- sequence_data$new("my_dataset")
    #'   metadata <- readr::read_tsv(rdataset_example("sample-metadata.tsv"),
    #'    col_names = TRUE, show_col_types = FALSE)
    #'   dataset$add_metadata(metadata)
    #'   dataset$get_metadata()
    #'
    #' @return data.frame()
    get_metadata = function() {
      private$metadata
    },

    #' @description
    #' Get the number of bins in the dataset
    #' @param type a string indicating the type of clusters. Default = "otu".
    #' @examples
    #'
    #'  dataset <- sequence_data$new("my_dataset")
    #'  bin_ids <- c("bin1", "bin2", "bin3")
    #'  abundances <- c(110, 525, 80)
    #'  dataset$assign_bins(bin_ids, abundances)
    #'  dataset$get_num_bins()
    #'
    #' @return An integer
    get_num_bins = function(type = "otu") {
      get_num_bins(self$data, type)
    },

    #' @description
    #' Get number of samples in the dataset
    #' @examples
    #'
    #' dataset <- sequence_data$new("my_dataset")
    #' bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
    #' samples <- c("sample1", "sample2", "sample5", "sample1",
    #'  "sample3", "sample1")
    #' sample_abundances <- c(10, 100, 1, 500, 25, 80)
    #' dataset$assign_bins(bin_ids, sample_abundances, samples)
    #' dataset$get_num_samples()
    #'
    #' @return An integer
    get_num_samples = function() {
      get_num_samples(self$data)
    },

    #' @description
    #' Get the number of sequences in the dataset
    #' @param distinct Boolean. When distinct is TRUE the number of unique
    #' sequence is returned.
    #' @param sample The name of the sample you want number of sequences for,
    #'  optional
    #' @return An integer
    get_num_sequences = function(distinct = FALSE, sample = NULL) {
      if (is.null(sample)) {
        sample <- ""
      }

      get_num_sequences(self$data, distinct, sample)
    },

    #' @description
    #' Get the number of treatments in the dataset
    #' @return An integer
    get_num_treatments = function() {
      get_num_treatments(self$data)
    },

    #' @description
    #' Get data.frame containing the oligo data
    #' tag, oligo, diffs, oligo(optional), diffs(optional), sample(optional)
    #' @return data.frame
    get_oligos = function() {
      # TODO
    },

    #' @description
    #' Get data frame containing sequence bin assignments
    #' @param type a string indicating the type of clusters. Options
    #' include: "otu", "asv", or "phylotype". Default = "otu".
    #' @examples
    #'   dataset <- sequence_data$new("my_dataset")
    #'   bin_ids <- c("bin1", "bin2", "bin3")
    #'   bin_abundances <- c(111, 525, 80)
    #'   dataset$assign_bins(bin_ids, bin_abundances)
    #'
    #'   # (rabund) bins would look like:
    #'   # bin1  bin2  bin3
    #'   # 110   525   80
    #'
    #'   rabund <- dataset$get_rabund()
    #'
    #'   #  rabund$bin_id  rabund$abundance
    #'   #  bin1           111
    #'   #  bin2           525
    #'   #  bin3           80
    #'
    #' @return data.frame
    get_rabund = function(type = "otu") {
      get_rabund(self$data, type)
    },

    #' @description
    #' Get data.frame containing information about the references used in the
    #' analysis of the dataset
    #' @examples
    #'   dataset <- sequence_data$new("my_dataset")
    #'   reference <- readr::read_csv(rdataset_example("references.csv"))
    #'   dataset$add_references(reference = reference)
    #'   dataset$get_references()
    #'
    #' @return data.frame()
    get_references = function() {
      private$references
    },

    #' @description
    #' Get names of samples in the dataset
    #' @return A character vector
    get_samples = function() {
      get_samples(self$data)
    },

    #' @description
    #' Get a summary of the samples and treatments in the dataset
    #' @param silent Default = FALSE, meaning print sample summary
    #' @return list
    get_sample_summary = function(silent = FALSE) {
      if (get_num_samples(self$data) != 0) {
        sample_totals <- get_sample_totals(self$data)
        sample_names <- self$get_samples()

        if (!silent) {
          cat("Sample   Total:\n")
          for (i in seq_along(sample_names)) {
            cat(paste(sample_names[i], sample_totals[i], sep = "\t"), "\n")
          }

          if (self$get_num_treatments() != 0) {
            treatment_names <- self$get_treatments()
            treatment_totals <- get_treatment_totals(self$data)
            cat("\n\n")
            cat("Treatment   Total:\n")
            for (i in seq_along(treatment_names)) {
              cat(
                paste(treatment_names[i], treatment_totals[i], sep = "\t"),
                "\n"
              )
            }
          }
        }

        if (self$get_num_treatments() != 0) {
          treatment_names <- self$get_treatments()
          treatment_totals <- get_treatment_totals(self$data)
          return(list(
            data.frame(sample = sample_names, total = sample_totals),
            data.frame(treatment = treatment_names, total = treatment_totals)
          ))
        } else {
          return(list(data.frame(sample = sample_names, total = sample_totals)))
        }
      } else {
        cli::cli_alert("Your dataset does not include sample data, ignoring.")
      }
      list()
    },

    #' @description
    #' Get phylo tree relating the samples in your dataset.
    #' @examples
    #'
    #'  dataset <- sequence_data$new("my_dataset")
    #'  df <- read_mothur_shared(rdataset_example("final.opti_mcc.shared"))
    #'  dataset$assign_bins(df$bin_id, df$abundance, df$sample)
    #'  tree <- ape::read.tree(rdataset_example(
    #'  "final.opti_mcc.jclass.ave.tre"))
    #'  dataset$add_sample_tree(tree)
    #'  dataset$get_sample_tree()
    #'
    get_sample_tree = function() {
      if (!is.null(private$sample_tree)) {
        # prune tree if needed
        # samples in tree and not in dataset
        extra_samples <- setdiff(
          private$sample_tree$tip.label,
          self$get_samples()
        )

        if (length(extra_samples) != 0) {
          # if tree contains "extra" samples, prune the tree
          private$sample_tree <- drop.tip(private$sample_tree,
            tip = extra_samples
          )
        }
      }
      private$sample_tree
    },

    #' @description
    #' Get report containing the scrapped sequences / bins -
    #' ids, trash_codes
    #' @return list of data.frames
    get_scrap_report = function() {
      results <- list()
      list_names <- c("sequence_scrap_report")
      scrap_sequence_report <- get_scrap_report(self$data, "sequence")
      results[[1]] <- scrap_sequence_report
      if (get_num_bins(self$data, "otu") != 0) {
        results[[2]] <- get_scrap_report(self$data, "otu")
        list_names <- c(list_names, "otu_scrap_report")
      }
      if (get_num_bins(self$data, "asv") != 0) {
        results[[3]] <- get_scrap_report(self$data, "asv")
        list_names <- c(list_names, "asv_scrap_report")
      }
      if (get_num_bins(self$data, "phylotype") != 0) {
        results[[4]] <- get_scrap_report(self$data, "phylotype")
        list_names <- c(list_names, "phylotype_scrap_report")
      }
      names(results) <- list_names
      results
    },

    #' @description
    #' Get the names of the sequences in the dataset
    #' @param sample a string containing the name of the sample you
    #' would like sequence names for. Default = NULL, meaning all samples in
    #' the dataset.
    #' @examples
    #'
    #' sequence_names <- c("seq1", "seq1", "seq1", "seq2", "seq2",
    #'              "seq2", "seq3", "seq3", "seq4")
    #' samples <- c("sample2", "sample3", "sample4",
    #'             "sample2", "sample3", "sample4",
    #'            "sample2", "sample3",
    #'            "sample4")
    #' abundances <- c(250, 400, 500,
    #'                25, 40, 50,
    #'                25, 25,
    #'                4)
    #'
    #' dataset <- sequence_data$new("my_dataset")
    #' dataset$assign_sequence_abundance(data = NULL,
    #'                                 names, abundances, samples)
    #' dataset$get_sequence_names()
    #'
    #' @return vector of sequence names
    get_sequence_names = function(sample = NULL) {
      if (is.null(sample)) {
        sample <- ""
      }

      get_sequence_names(self$data, sample)
    },

    #' @description
    #' Get data.frame sequence report data. Sequence report data includes: start
    #' positions, end positions, number of bases, number of ambiguous bases,
    #' length of longest homopolymer, and the number of N's.
    #' @return data.frame
    get_sequence_report = function() {
      get_sequence_report(self$data)
    },

    #' @description
    #' Get summary of the sequence reports
    #' @param silent Default = FALSE, meaning print summaries
    #' @return list of data.frames
    get_sequence_summary = function(silent = FALSE) {
      results <- get_sequence_summary(self$data)

      if (nrow(private$contigs_data) != 0) {
        results[["contigs_summary"]] <- private$summarize(
          self$get_contigs_assembly_report()
        )

        # if you have summary results to print
        if (!silent) {
          print(results[["contigs_summary"]])
          cat("Unique seqs:\t", self$get_num_sequences(TRUE), "\n")
          cat("Total seqs:\t", self$get_num_sequences(), "\n")
        }
      }

      # if you have alignment data, then print
      if (nrow(private$alignment_data) != 0) {
        results[["alignment_summary"]] <- private$summarize(
          self$get_alignment_report()
        )

        # if you have summary results to print
        if (!silent) {
          print(results[["alignment_summary"]])
          cat("Unique seqs:\t", self$get_num_sequences(TRUE), "\n")
          cat("Total seqs:\t", self$get_num_sequences(), "\n")
        }
      }

      # if you have summary results to print
      if (!all(self$get_sequences() == "") && (!silent)) {
        # if you have summary results to print
        if (!silent) {
          print(results[[1]])
          cat("Unique seqs:\t", self$get_num_sequences(TRUE), "\n")
          cat("Total seqs:\t", self$get_num_sequences(), "\n")
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

      return(results)
    },

    #' @description
    #' Get report containing the sequence taxonomy table -
    #' ids, taxonomy by levels
    #' @examples
    #' names <- c("seq1", "seq2", "seq3", "seq4")
    #' tax1 <- "Bacteria(100);Bacteroidetes(95);Bacteroidia(90);"
    #' tax2 <- "Bacteria(100);Proteobacteria(89);Betaproteobacteria(85);"
    #' tax3 <- "Bacteria(100);Firmicutes(99);Bacilli(90);"
    #' tax4 <- "Bacteria(100);Proteobacteria(87);Gammaproteobacteria(82);"
    #' taxonomies <- c(tax1, tax2, tax3, tax4)
    #'
    #' dataset <- sequence_data$new("my_dataset")
    #' dataset$assign_sequence_taxonomy(names, taxonomies)
    #' dataset$get_sequence_taxonomy_report()
    #'
    #' @return data.frame
    get_sequence_taxonomy_report = function() {
      get_sequence_taxonomy_report(self$data)
    },

    #' @description
    #' Get phylo tree relating the sequences in your dataset.
    #' @examples
    #'
    #'  dataset <- sequence_data$new("my_dataset")
    #'  tree <- ape::read.tree(rdataset_example("final.phylip.tre"))
    #'  dataset$add_sequence_tree(tree)
    #'  dataset$get_sequence_tree()
    #'
    get_sequence_tree = function() {
      if (!is.null(private$sequence_tree)) {
        # prune tree if needed
        # seqs in tree and not in dataset
        extra_seqs <- setdiff(
          private$sequence_tree$tip.label,
          self$get_sequence_names()
        )

        if (length(extra_seqs) != 0) {
          # if tree contains "extra" names, prune the tree
          private$sequence_tree <- drop.tip(private$sequence_tree,
            tip = extra_seqs
          )
        }
      }
      private$sequence_tree
    },

    #' @description
    #' Get vector containing sequence nucleotide strings
    #' @param sample String, name of sample
    get_sequences = function(sample = NULL) {
      if (is.null(sample)) {
        sample <- ""
      }

      get_sequences(self$data, sample)
    },

    #' @description
    #' Get names of treatments in the dataset
    #' @return A character vector
    get_treatments = function() {
      get_treatments(self$data)
    },

    #' @description
    #' Determine if a sample is present in the dataset
    #' @param sample String, Name of sample
    #' @return Boolean
    has_sample = function(sample) {
      has_sample(self$data, sample)
    },

    #' @description
    #' Determine if the dataset is aligned
    #' @return bool
    is_aligned = function() {
      is_aligned(self$data)
    },

    #' @description
    #' Remove contaminants from the dataset
    #' @param contaminants vector of strings containing the taxonomies you would
    #' like to remove
    #' @examples
    #' dataset <- read_mothur(fasta = rdataset_example("final.fasta"),
    #'                       count = rdataset_example("final.count_table"),
    #'                       taxonomy = rdataset_example("final.taxonomy"),
    #'                       design = rdataset_example("mouse.time.design"),
    #'                       otu_list = rdataset_example("final.opti_mcc.list"),
    #'                       dataset_name = "miseq_sop")
    #'
    #' contaminants <- c("Chloroplast", "Mitochondria", "unknown", "Archaea",
    #'  "Eukaryota")
    #'
    #' dataset$remove_lineages(contaminants)
    #'
    remove_lineages = function(contaminants) {
      remove_lineages(self$data, contaminants, "contaminant")
      invisible(self)
    },


    #' @description
    #' Remove samples from the dataset
    #' @param samples vector of strings containing the names of the samples to
    #' @examples
    #' dataset <- read_mothur(fasta = rdataset_example("final.fasta"),
    #'                       count = rdataset_example("final.count_table"),
    #'                       taxonomy = rdataset_example("final.taxonomy"),
    #'                       design = rdataset_example("mouse.time.design"),
    #'                       otu_list = rdataset_example("final.opti_mcc.list"),
    #'                       dataset_name = "miseq_sop")
    #'
    #' dataset$get_num_samples()
    #'
    #' # To remove samples 'F3D0' and 'F3D1'
    #'
    #' dataset$remove_samples(c("F3D0", "F3D1"))
    #'
    #' dataset$get_num_samples()
    #'
    remove_samples = function(samples) {
      remove_samples(self$data, samples)
      invisible(self)
    }
  ),
  private = list(
    metadata = data.frame(),
    references = data.frame(),
    alignment_data = data.frame(),
    contigs_data = data.frame(),
    sequence_tree = NULL,
    sample_tree = NULL,
    processors = 1,
    finalize = function() {},

    # summarize numeric columns in a dataframe
    summarize = function(report) {
      # rcpp function - calls summary.cpp
      report_summary <- summarize_reports(
        report[sapply(report, is.numeric)],
        get_sequence_abundances(self$data),
        private$processors
      )
    },
    fill_required_param = function(param, data, default_value) {
      data_names <- names(data)

      # required
      if (is.null(param)) {
        param <- default_value
      }

      if (length(param) == 1) {
        if (param %in% data_names) {
          param <- data[[param]]
        } else {
          abort_missing_column(param)
        }
      }

      param
    },
    fill_optional_param = function(param, data, default_value) {
      data_names <- names(data)

      # optional
      if (!is.null(param)) {
        if (length(param) == 1) {
          if (param %in% data_names) {
            param <- data[[param]]
          } else {
            abort_missing_column(param)
          }
        }
      } else {
        # look for default column
        param <- default_value
        if (param %in% data_names) {
          param <- data[[param]]
        } else {
          param <- ""
        }
      }

      param
    }
  )
)
