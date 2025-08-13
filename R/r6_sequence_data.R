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

    #' @field data pointer to 'Dataset' (Rcpp Module). This allows package
    #' developers an easy access point to the underlying C++ code with
    #' additional functionality.
    data = NULL,

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
        self$data <- new(Dataset, name, processors)
      } else {
        # copy of dataset backend
        self$data <- new(Dataset, dataset$data)
        self$data$processors <- processors
        # assign new name
        if (name != "") {
          self$data$dataset_name <- name
        }

        # copy metadata
        private$metadata <- dataset$get_metadata()
      }

      invisible(self)
    },

    #' @description
    #' Get summary of sequence data
    print = function() {
      if (self$get_dataset_name() != "") {
        cat(self$get_dataset_name())
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
        abort_incorrect_type("data.frame", class(metadata)[1])
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
    #' reference <- readr::read_csv(rdataset_example("references.csv"))
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
          abort_incorrect_type("data.frame", class(reference)[1])
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
    #' Add new sequence data
    #' @param names a vector of sequence names
    #' @param sequences a vector of sequence data
    #' @param comments a vector of sequence comments, (optional)
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
    #'   sequences <- read_fasta(rdataset_example("final.fasta"))
    #'   dataset$add_sequences(sequences$names, sequences$sequences)
    #'
    #' # With the additional parameters to add information about the reference
    #' # You can also add references using the 'add_references' function.
    #'
    #' url <- "https://mothur.org/wiki/silva_reference_files/"
    #'
    #' dataset <- sequence_data$new("my_dataset")
    #' dataset$add_sequences(sequences$names, sequences$sequences,
    #' reference_name = "silva.bacteria.fasta",
    #' reference_note = "alignment by mothur2 v1.0 using default options",
    #' reference_version = "1.38.1", reference_url = url)
    #'
    add_sequences = function(names, sequences = NULL, comments = NULL,
                             reference_name = NULL,
                             reference_version = NULL,
                             reference_note = NULL,
                             reference_url = NULL) {
      if (is.null(comments)) {
        comments <- rep("", length(names))
      }
      if (is.null(sequences)) {
        sequences <- rep("", length(names))
      }
      self$data$add_sequences(names, sequences, comments)

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
    #' Add bin data
    #' @param bin_ids a vector of bin labels
    #' @param abundances a vector of abundances (optional). You must provide
    #'  either abundances or seq_ids.
    #' @param samples a vector of sample assignments (optional)
    #' @param seq_ids a vector of sequence names (optional) You must provide
    #'  either abundances or seq_ids.
    #' @param type a string indicating the type of bin assignments. Options
    #' include: "otu", "asv", or "phylotype". Default = "otu".
    #' @examples
    #'
    #'   # bin_ids  seq_ids
    #'   #  bin1     seq1
    #'   #  bin1     seq2
    #'   #  bin1     seq4
    #'   #  bin2     seq3
    #'   #  bin2     seq6
    #'   #  bin3     seq5
    #'
    #'   # To assign sequences to bins:
    #'
    #'   dataset <- sequence_data$new("my_dataset")
    #'   seq_ids <- c("seq1", "seq2", "seq4", "seq3", "seq6", "seq5")
    #'   bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
    #'   dataset$assign_bins(bin_ids, seq_ids = seq_ids)
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
    #'   dataset$assign_bins(bin_ids, abundances)
    #'
    #'   # bins would look like:
    #'   #            bin1             bin2        bin3
    #'   # (rabund)   110              525         80
    #'
    #'   # To add abundance bin assignments parsed by sample:
    #'
    #'   # bin_ids  samples  abundances
    #'   #  bin1     sample1        10
    #'   #  bin1     sample2        100
    #'   #  bin1     sample5        1
    #'   #  bin2     sample1        500
    #'   #  bin2     sample3        25
    #'   #  bin3     sample1        80
    #'
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
    #'   dataset$assign_bins(bin_ids, abundances, samples, seq_ids)
    #'
    #'   # bins would look like:
    #'   #            bin1             bin2        bin3
    #'   # (list)     seq1,seq2,seq4   seq3,seq6   seq5
    #'   # (rabund)   716              85          65
    #'
    #'   dataset$get_shared()
    #'
    assign_bins = function(bin_ids, abundances = NULL,
                           samples = NULL, seq_ids = NULL, type = "otu") {
      if (is.null(abundances) && is.null(seq_ids)) {
        cli::cli_abort("[ERROR]: You must provide either
                           abundances or seq_ids.")
      }

      if (is.null(samples)) {
        samples <- rep("", length(bin_ids))
      }
      if (is.null(abundances)) {
        abundances <- rep(0, length(bin_ids))
      }
      if (is.null(seq_ids)) {
        seq_ids <- rep("", length(bin_ids))
      }

      self$data$assign_bins(
        bin_ids, abundances,
        samples, seq_ids, type
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
    #' @param type a string indicating the type of clusters. Options
    #' include: "otu", "asv", or "phylotype". Default = "otu".
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
      if (self$data$get_num_bins(type) == 0) {
        cli::cli_abort("[ERROR]: No bin data for type " + type + ", please
                          assign bins using the 'assign_bins' function then
                       try again.")
      }
      self$data$assign_bin_taxonomy(bin_ids, taxonomies, type)

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
    #' @param names a vector of sequence names
    #' @param abundances a vector of sequence abundances
    #' @param samples a vector of sample assignments (optional)
    #' @param treatments a vector of treatment assignments (optional)
    #' @examples
    #'
    #'  # mothur count file
    #'  # Representative_Sequence     total   sample2	sample3	sample4
    #'  # seq1	1150	250	400	500
    #'  # seq2	115	25	40	50
    #'  # seq3	50	25	25	0
    #'  # seq4	4	0	0	4
    #'
    #' # inputted as a sample table
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
    #'
    #' dataset <- sequence_data$new("my_dataset")
    #' unique_names <- unique(names)
    #' dataset$add_sequences(unique_names)
    #' dataset$assign_sequence_abundance(names, abundances, samples)
    #'
    assign_sequence_abundance = function(names, abundances, samples = NULL,
                                         treatments = NULL) {
      if (length(names) != length(abundances)) {
        cli::cli_abort("[ERROR]: The names and abundances must be the same
                         length.")
      }

      # sanity check, make sure names are present in dataset
      unique_names <- sort(unique(names))
      dataset_names <- sort(self$get_ids())

      if (!identical(unique_names, dataset_names)) {
        cli::cli_abort("[ERROR]: You must provide assignments for all
                          sequences in your dataset.")
      }
      unique_names <- c()
      dataset_names <- c()

      if (is.null(samples)) {
        samples <- rep("", length(names))
      }

      if (is.null(treatments)) {
        treatments <- rep("", length(names))
      }

      self$data$assign_sequence_abundance(
        names, abundances,
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
      self$data$assign_sequence_taxonomy(names, taxonomies)

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
    #'
    #' dataset <- sequence_data$new("my_dataset")
    #' unique_names <- unique(names)
    #' dataset$add_sequences(unique_names)
    #' dataset$assign_sequence_abundance(names, abundances, samples)
    #'
    #' treatments <- c("early", "early", "late")
    #' dataset$assign_treatments(unique(samples), treatments)
    #'
    assign_treatments = function(samples, treatments) {
      self$data$assign_treatments(samples, treatments)
      invisible(self)
    },

    #' @description
    #' Remove all data from dataset
    clear = function() {
      self$data$clear()
      metadata <- data.frame()
      references <- data.frame()
      invisible(self)
    },

    #' @description
    #' Get List containing dataset
    #' @return List
    export = function() {
      self$data$export()
    },

    #' @description
    #' Get report containing align report table -
    #' ids, search_scores, sim_scores, longest_insert
    #' @return data.frame
    get_align_report = function() {
      if (self$data$has_align_data) {
        return(self$data$get_align_report())
      }
      data.frame()
    },

    #' @description
    #' Get report containing contigs report table -
    #' ids, lengths, ostarts, oends, olengths, mismatches, numns, ee
    #' @return data.frame
    get_contigs_report = function() {
      if (self$data$has_contigs_data) {
        return(self$data$get_contigs_report())
      }
      data.frame()
    },

    #' @description
    #' Get count table returns data.frame containing:
    #' ids, abundances, sample(optional), treatment(optional)
    #' This table represents mothur's count and design files.
    #' @return data.frame
    get_count_table = function() {
      self$data$get_sequence_abundance_table()
    },

    #' @description
    #' Get dataset name
    #' @return String
    get_dataset_name = function() {
      self$data$dataset_name
    },

    #' @description
    #' Get names of sequences in the dataset
    #' @param sample String, name of sample
    #' @return vector of ids
    get_ids = function(sample = NULL) {
      if (is.null(sample)) {
        sample <- ""
      }

      self$data$get_names(sample)
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
    #'   dataset$assign_bins(bin_ids, sequence_abundances, seq_ids = seq_ids)
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
      self$data$get_list(type)
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
    #' @param type a string indicating the type of clusters. Options
    #' include: "otu", "asv", or "phylotype". Default = "otu".
    #' @return An integer
    get_num_bins = function(type = "otu") {
      self$data$get_num_bins(type)
    },

    #' @description
    #' Get number of samples in the dataset
    #' @return An integer
    get_num_samples = function() {
      self$data$get_num_samples()
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

      if (distinct) {
        self$data$get_unique_total(sample)
      } else {
        self$data$get_total(sample)
      }
    },

    #' @description
    #' Get the number of treatments in the dataset
    #' @return An integer
    get_num_treatments = function() {
      self$data$get_num_treatments()
    },

    #' @description
    #' Get data.frame containing the oligo data
    #' tag, oligo, diffs, oligo(optional), diffs(optional), sample(optional)
    #' @return data.frame
    get_oligos = function() {
      # TODO
    },

    #' @description
    #' Get report containing the bin taxonomy table -
    #' ids, taxonomy by levels
    #' @param type a string indicating the type of clusters. Options
    #' include: "otu", "asv", or "phylotype". Default = "otu".
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
      self$data$get_bin_taxonomy_report(type)
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
      self$data$get_rabund(type)
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
      self$data$get_samples()
    },

    #' @description
    #' Get summary of count table. sample, sample_totals
    #' @param silent Default = FALSE, meaning print sample summary
    #' @return list
    get_sample_summary = function(silent = FALSE) {
      if (self$data$get_num_samples() != 0) {
        sample_totals <- self$data$get_sample_totals()
        sample_names <- self$get_samples()

        if (!silent) {
          cat("Sample   Total:\n")
          for (i in seq_along(sample_names)) {
            cat(paste(sample_names[i], sample_totals[i], sep = "\t"), "\n")
          }

          if (self$get_num_treatments() != 0) {
            treatment_names <- self$get_treatments()
            treatment_totals <- self$data$get_treatment_totals()
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
          treatment_totals <- self$data$get_treatment_totals()
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
    #' Get report containing the scrapped sequences / bins -
    #' ids, trash_codes
    #' @return list of data.frames
    get_scrap_report = function() {
      results <- list()
      list_names <- c("sequence_scrap_report")
      scrap_sequence_report <- self$data$get_scrap_report("sequence")
      results[[1]] <- scrap_sequence_report
      if (self$data$get_num_bins("otu") != 0) {
        results[[2]] <- self$data$get_scrap_report("otu")
        list_names <- c(list_names, "otu_scrap_report")
      }
      if (self$data$get_num_bins("asv") != 0) {
        results[[3]] <- self$data$get_scrap_report("asv")
        list_names <- c(list_names, "asv_scrap_report")
      }
      if (self$data$get_num_bins("phylotype") != 0) {
        results[[4]] <- self$data$get_scrap_report("phylotype")
        list_names <- c(list_names, "phylotype_scrap_report")
      }
      names(results) <- list_names
      results
    },

    #' @description
    #' Get data.frame sequence report data. Sequence report data includes: start
    #' positions, end positions, number of bases, number of ambiguous bases,
    #' length of longest homopolymer, and the number of N's.
    #' @return data.frame
    get_sequence_report = function() {
      self$data$get_sequence_report()
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
      self$data$get_sequence_taxonomy_report()
    },

    #' @description
    #' Get summary of the sequence reports
    #' @param silent Default = FALSE, meaning print summaries
    #' @return list of data.frames
    get_sequence_summary = function(silent = FALSE) {
      results <- self$data$get_sequence_summary()

      if (length(results) == 0) {
        cli::cli_alert("Your dataset does not include sequence data, ignoring.")
        return(list())
      }

      results_row_names <- c(
        "Minimum:", "2.5%-tile:", "25%-tile:",
        "Median:   ", "75%-tile:", "97.5%-tile:",
        "Maximum:", "Mean:      "
      )

      rownames(results[[1]]) <- results_row_names

      if (self$data$has_contigs_data) {
        rownames(results$contigs_summary) <- results_row_names

        # if you have summary results to print
        if ((length(results$contigs_summary$ostarts) == 8) && (!silent)) {
          cat("\t\toverlap_start\toverlap_end\tlength\toverlap_length\t")
          cat("mismatches\tnum_ns\tnumseqs\n")
          cat(paste(
            results_row_names[1], results$contigs_summary$ostarts[1],
            results$contigs_summary$oends[1],
            results$contigs_summary$lengths[1],
            results$contigs_summary$olengths[1],
            results$contigs_summary$mismatches[1],
            results$contigs_summary$numns[1],
            results$contigs_summary$numseqs[1],
            sep = "\t"
          ), "\n")
          cat(paste(
            results_row_names[2], results$contigs_summary$ostarts[2],
            results$contigs_summary$oends[2],
            results$contigs_summary$lengths[2],
            results$contigs_summary$olengths[2],
            results$contigs_summary$mismatches[2],
            results$contigs_summary$numns[2],
            results$contigs_summary$numseqs[2],
            sep = "\t"
          ), "\n")
          cat(paste(
            results_row_names[3], results$contigs_summary$ostarts[3],
            results$contigs_summary$oends[3],
            results$contigs_summary$lengths[3],
            results$contigs_summary$olengths[3],
            results$contigs_summary$mismatches[3],
            results$contigs_summary$numns[3],
            results$contigs_summary$numseqs[3],
            sep = "\t"
          ), "\n")
          cat(paste(
            results_row_names[4], results$contigs_summary$ostarts[4],
            results$contigs_summary$oends[4],
            results$contigs_summary$lengths[4],
            results$contigs_summary$olengths[4],
            results$contigs_summary$mismatches[4],
            results$contigs_summary$numns[4],
            results$contigs_summary$numseqs[4],
            sep = "\t"
          ), "\n")
          cat(paste(
            results_row_names[5], results$contigs_summary$ostarts[5],
            results$contigs_summary$oends[5],
            results$contigs_summary$lengths[5],
            results$contigs_summary$olengths[5],
            results$contigs_summary$mismatches[5],
            results$contigs_summary$numns[5],
            results$contigs_summary$numseqs[5],
            sep = "\t"
          ), "\n")
          cat(paste(
            results_row_names[6], results$contigs_summary$ostarts[6],
            results$contigs_summary$oends[6],
            results$contigs_summary$lengths[6],
            results$contigs_summary$olengths[6],
            results$contigs_summary$mismatches[6],
            results$contigs_summary$numns[6],
            results$contigs_summary$numseqs[6],
            sep = "\t"
          ), "\n")
          cat(paste(
            results_row_names[7], results$contigs_summary$ostarts[7],
            results$contigs_summary$oends[7],
            results$contigs_summary$lengths[7],
            results$contigs_summary$olengths[7],
            results$contigs_summary$mismatches[7],
            results$contigs_summary$numns[7],
            results$contigs_summary$numseqs[7],
            sep = "\t"
          ), "\n")
          cat(paste(
            results_row_names[8], results$contigs_summary$ostarts[8],
            results$contigs_summary$oends[8],
            results$contigs_summary$lengths[8],
            results$contigs_summary$olengths[8],
            results$contigs_summary$mismatches[8],
            results$contigs_summary$numns[8],
            sep = "\t"
          ), "\n")
          cat("Unique seqs:\t", self$get_num_sequences(TRUE), "\n")
          cat("Total seqs:\t", self$get_num_sequences(), "\n")
        }
      }

      # if you have alignment data, then print
      if (self$data$has_align_data) {
        rownames(results$align_summary) <- results_row_names

        if ((length(results$align_summary$search_scores) == 8) && (!silent)) {
          cat("\t\tsearch_scores\tsim_scores\tlongest_inserts\n")
          cat(paste(
            results_row_names[1], results$align_summary$search_scores[1],
            results$align_summary$sim_scores[1],
            results$align_summary$longest_inserts[1],
            sep = "\t"
          ), "\n")
          cat(paste(
            results_row_names[2], results$align_summary$search_scores[2],
            results$align_summary$sim_scores[2],
            results$align_summary$longest_inserts[2],
            sep = "\t"
          ), "\n")
          cat(paste(
            results_row_names[3], results$align_summary$search_scores[3],
            results$align_summary$sim_scores[3],
            results$align_summary$longest_inserts[3],
            sep = "\t"
          ), "\n")
          cat(paste(
            results_row_names[4], results$align_summary$search_scores[4],
            results$align_summary$sim_scores[4],
            results$align_summary$longest_inserts[4],
            sep = "\t"
          ), "\n")
          cat(paste(
            results_row_names[5], results$align_summary$search_scores[5],
            results$align_summary$sim_scores[5],
            results$align_summary$longest_inserts[5],
            sep = "\t"
          ), "\n")
          cat(paste(
            results_row_names[6], results$align_summary$search_scores[6],
            results$align_summary$sim_scores[6],
            results$align_summary$longest_inserts[6],
            sep = "\t"
          ), "\n")
          cat(paste(
            results_row_names[7], results$align_summary$search_scores[7],
            results$align_summary$sim_scores[7],
            results$align_summary$longest_inserts[7],
            sep = "\t"
          ), "\n")
          cat(paste(
            results_row_names[8], results$align_summary$search_scores[8],
            results$align_summary$sim_scores[8],
            results$align_summary$longest_inserts[8],
            sep = "\t"
          ), "\n")
          cat("Unique seqs:\t", self$get_num_sequences(TRUE), "\n")
          cat("Total seqs:\t", self$get_num_sequences(), "\n")
        }
      }

      # if you have summary results to print
      if (self$data$has_sequence_strings() && (!silent)) {
        cat("\t\tstart\tend\tlength\tambigs\tpolymer\tnum_ns\tnumseqs\n")
        cat(
          paste(results_row_names[1], results[[1]]$starts[1],
            results[[1]]$ends[1],
            results[[1]]$nbases[1], results[[1]]$ambigs[1],
            results[[1]]$polymers[1], results[[1]]$numns[1],
            results[[1]]$numseqs[1],
            sep = "\t"
          ),
          "\n"
        )
        cat(
          paste(results_row_names[2], results[[1]]$starts[2],
            results[[1]]$ends[2],
            results[[1]]$nbases[2], results[[1]]$ambigs[2],
            results[[1]]$polymers[2], results[[1]]$numns[2],
            results[[1]]$numseqs[2],
            sep = "\t"
          ),
          "\n"
        )
        cat(
          paste(results_row_names[3], results[[1]]$starts[3],
            results[[1]]$ends[3],
            results[[1]]$nbases[3], results[[1]]$ambigs[3],
            results[[1]]$polymers[3], results[[1]]$numns[3],
            results[[1]]$numseqs[3],
            sep = "\t"
          ),
          "\n"
        )
        cat(
          paste(results_row_names[4], results[[1]]$starts[4],
            results[[1]]$ends[4],
            results[[1]]$nbases[4], results[[1]]$ambigs[4],
            results[[1]]$polymers[4], results[[1]]$numns[4],
            results[[1]]$numseqs[4],
            sep = "\t"
          ),
          "\n"
        )
        cat(
          paste(results_row_names[5], results[[1]]$starts[5],
            results[[1]]$ends[5],
            results[[1]]$nbases[5], results[[1]]$ambigs[5],
            results[[1]]$polymers[5], results[[1]]$numns[5],
            results[[1]]$numseqs[5],
            sep = "\t"
          ),
          "\n"
        )
        cat(
          paste(results_row_names[6], results[[1]]$starts[6],
            results[[1]]$ends[6],
            results[[1]]$nbases[6], results[[1]]$ambigs[6],
            results[[1]]$polymers[6], results[[1]]$numns[6],
            results[[1]]$numseqs[6],
            sep = "\t"
          ),
          "\n"
        )
        cat(
          paste(results_row_names[7], results[[1]]$starts[7],
            results[[1]]$ends[7],
            results[[1]]$nbases[7], results[[1]]$ambigs[7],
            results[[1]]$polymers[7], results[[1]]$numns[7],
            results[[1]]$numseqs[7],
            sep = "\t"
          ),
          "\n"
        )
        cat(paste(
          results_row_names[8], results[[1]]$starts[8], results[[1]]$ends[8],
          results[[1]]$nbases[8], results[[1]]$ambigs[8],
          results[[1]]$polymers[8], results[[1]]$numns[8],
          sep = "\t"
        ), "\n")
        cat("Unique seqs:\t", self$get_num_sequences(TRUE), "\n")
        cat("Total seqs:\t", self$get_num_sequences(), "\n")
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
    #' Get vector containing sequence nucleotide strings
    #' @param sample String, name of sample
    get_sequences = function(sample = NULL) {
      if (is.null(sample)) {
        sample <- ""
      }

      self$data$get_sequences(sample)
    },

    #' @description
    #' Get data frame containing sequence bin assignments by sample
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
    #'   shared <- dataset$get_shared()
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
    get_shared = function(type = "otu") {
      self$data$get_shared(type)
    },

    #' @description
    #' Get names of treatments in the dataset
    #' @return A character vector
    get_treatments = function() {
      self$data$get_treatments()
    },

    #' @description
    #' Determine if a sample is present in the dataset
    #' @param sample String, Name of sample
    #' @return Boolean
    has_sample = function(sample) {
      self$data$has_sample(sample)
    },

    #' @description
    #' Determine if the dataset is aligned
    #' @return bool
    is_aligned = function() {
      self$data$is_aligned
    },

    #' @description
    #' Remove contaminants from the dataset
    #' @param contaminants vector of strings containing the taxonomies you would
    #' like to remove
    #' @examples
    #' dataset <- read_mothur(fasta = rdataset_example("final.fasta"),
    #'                        count = rdataset_example("final.count_table"),
    #'                        taxonomy = rdataset_example("final.taxonomy"),
    #'                        design = rdataset_example("mouse.time.design"),
    #'                        list = rdataset_example("final.opti_mcc.list"),
    #'                        dataset_name = "miseq_sop")
    #'
    #' contaminants <- c("Chloroplast", "Mitochondria", "unknown", "Archaea",
    #'  "Eukaryota")
    #'
    #' dataset$remove_lineages(contaminants)
    #'
    remove_lineages = function(contaminants) {
      self$data$remove_lineages(contaminants, "contaminant")
    },


    #' @description
    #' Remove samples from the dataset
    #' @param samples vector of strings containing the names of the samples to
    #' @examples
    #' dataset <- read_mothur(fasta = rdataset_example("final.fasta"),
    #'                        count = rdataset_example("final.count_table"),
    #'                        taxonomy = rdataset_example("final.taxonomy"),
    #'                        design = rdataset_example("mouse.time.design"),
    #'                        list = rdataset_example("final.opti_mcc.list"),
    #'                        dataset_name = "miseq_sop")
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
      self$data$remove_samples(samples)
    }
  ),
  private = list(
    metadata = data.frame(),
    references = data.frame(),

    # Clear sequences from dataset
    finalize = function() {
      self$clear()
    }
  )
)
