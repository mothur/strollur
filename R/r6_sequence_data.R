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
    #' @param name String, name of dataset
    #' @param fasta String name of FASTA file
    #' @param processors Integer, number of cores to use.
    #'  Default = all available
    #' @examples
    #'
    #' # to create an empty dataset, run the following:
    #'
    #' dataset <- sequence_data$new("soil")
    #'
    #' # to create a dataset containing sequences from your fasta file,
    #' #  run the following:
    #'
    #' dataset <- sequence_data$new(name = "soil",
    #'                              fasta = rdataset_example("test.fasta"))
    #'
    #' @return A new `sequence_data` object.
    initialize = function(name, fasta = NULL,
                          processors = parallelly::availableCores()) {
      self$data <- new(Dataset, name, processors)

      if (!is.null(fasta)) {
        fasta_data <- read_fasta(fasta)
        self$add_sequences(fasta_data$names, fasta_data$sequences)
      }
      invisible(self)
    },

    #' @description
    #' Get summary of sequence data
    print = function() {
      cat(self$get_dataset_name())
      cat(":\n\n")
      self$get_sequence_summary()
      cat("\n\n")
      self$get_sample_summary()
      if (self$get_num_sequences(TRUE) != -1) {
        cat(
          paste("\nNumber of unique seqs:", self$get_num_sequences(TRUE)),
          "\n"
        )
      } else {
        cat("\n")
      }
      cat(
        paste("Total number of seqs:", self$get_num_sequences(), "\n"),
        "\n"
      )
    },

    #' @description
    #' Add new sequence data
    #' @param names a vector of sequence names
    #' @param sequences a vector of sequence data
    #' @param comments a vector of sequence comments, (optional)
    #' @examples
    #'
    #'   dataset <- sequence_data$new("my_dataset")
    #'   sequences <- read_fasta(rdataset_example("test.fasta"))
    #'   dataset$add_sequences(sequences$names, sequences$sequences)
    #'
    add_sequences = function(names, sequences = NULL, comments = NULL) {
      if (is.null(comments)) {
        comments <- rep("", length(names))
      }
      if (is.null(sequences)) {
        sequences <- rep("", length(names))
      }
      self$data$add_sequences(names, sequences, comments)

      invisible(self)
    },

    #' @description
    #' Add otu data
    #' @param otu_ids a vector of otu labels
    #' @param abundances a vector of abundances (optional). You must provide
    #'  either abundances or seq_ids.
    #' @param samples a vector of sample assignments (optional)
    #' @param seq_ids a vector of sequence names (optional) You must provide
    #'  either abundances or seq_ids.
    #' @examples
    #'
    #'   # otu_ids  seq_ids
    #'   #  otu1     seq1
    #'   #  otu1     seq2
    #'   #  otu1     seq4
    #'   #  otu2     seq3
    #'   #  otu2     seq6
    #'   #  otu3     seq5
    #'
    #'   # To assign sequences to otus:
    #'
    #'   dataset <- sequence_data$new("my_dataset")
    #'   seq_ids <- c("seq1", "seq2", "seq4", "seq3", "seq6", "seq5")
    #'   otu_ids <- c("otu1", "otu1", "otu1", "otu2", "otu2", "otu3")
    #'   dataset$assign_otus(otu_ids, seq_ids = seq_ids)
    #'
    #'   # otus would look like:
    #'   #            otu1             otu2        otu3
    #'   # (list)     seq1,seq2,seq4   seq3,seq6   seq5
    #'
    #'   # To add abundance only otu assignments:
    #'
    #'   dataset <- sequence_data$new("my_dataset")
    #'   otu_ids <- c("otu1", "otu2", "otu3")
    #'   abundances <- c(110, 525, 80)
    #'   dataset$assign_otus(otu_ids, abundances)
    #'
    #'   # otus would look like:
    #'   #            otu1             otu2        otu3
    #'   # (rabund)   110              525         80
    #'
    #'   # To add abundance otu assignments parsed by sample:
    #'
    #'   # otu_ids  samples  abundances
    #'   #  otu1     sample1        10
    #'   #  otu1     sample2        100
    #'   #  otu1     sample5        1
    #'   #  otu2     sample1        500
    #'   #  otu2     sample3        25
    #'   #  otu3     sample1        80
    #'
    #'   dataset <- sequence_data$new("my_dataset")
    #'   otu_ids <- c("otu1", "otu1", "otu1", "otu2", "otu2", "otu3")
    #'   samples <- c("sample1", "sample2", "sample5",
    #'    "sample1", "sample3", "sample1")
    #'   sample_abundances <- c(10, 100, 1, 500, 25, 80)
    #'   dataset$assign_otus(otu_ids, sample_abundances, samples)
    #'
    #'   # (shared) otus would look like:
    #'   # label  sample   otu1   otu2   otu3
    #'   # 0.03   sample1  10     500    80
    #'   # 0.03   sample2  100    0      0
    #'   # 0.03   sample3  0      25     0
    #'   # 0.03   sample5  1      0      0
    #'
    #'   # To assign sequences to otus with their abundances parsed by sample:
    #'
    #'   dataset <- sequence_data$new("my_dataset")
    #'   otu_ids <- c("otu1", "otu1", "otu1", "otu1", "otu1", "otu1",
    #'                "otu2", "otu2", "otu2",
    #'                "otu3", "otu3")
    #'   seq_ids <- c("seq1", "seq1", "seq1", "seq2", "seq4", "seq4",
    #'                "seq3", "seq3", "seq6",
    #'                "seq5", "seq5")
    #'   samples <- c("sample1", "sample2", "sample5",
    #'                "sample1", "sample3", "sample4",
    #'                "sample2", "sample3", "sample1",
    #'                "sample1", "sample6")
    #'   abundances <- c(10, 100, 1, 500, 25, 80, 20, 5, 60, 15, 50)
    #'   dataset$assign_otus(otu_ids, abundances, samples, seq_ids)
    #'
    #'   # otus would look like:
    #'   #            otu1             otu2        otu3
    #'   # (list)     seq1,seq2,seq4   seq3,seq6   seq5
    #'   # (rabund)   716              85          65
    #'
    #'   dataset$get_shared()
    #'
    assign_otus = function(otu_ids, abundances = NULL,
                           samples = NULL, seq_ids = NULL) {
      if (is.null(abundances) && is.null(seq_ids)) {
        cli::cli_abort("[ERROR]: You must provide either
                           abundances or seq_ids.")
      }

      if (is.null(samples)) {
        samples <- rep("", length(otu_ids))
      }
      if (is.null(abundances)) {
        abundances <- rep(0, length(otu_ids))
      }
      if (is.null(seq_ids)) {
        seq_ids <- rep("", length(otu_ids))
      }

      self$data$assign_otus(
        otu_ids, abundances,
        samples, seq_ids
      )

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
    #' @examples
    #'
    #' names <- c("seq1", "seq2", "seq3", "seq4")
    #' taxonomies <- c("Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;",
    #'               "Bacteria;Proteobacteria;Betaproteobacteria;Neisseriales;",
    #'                "Bacteria;Firmicutes;Bacilli;Lactobacillales;",
    #'            "Bacteria;Proteobacteria;Gammaproteobacteria;Pasteurellales;")
    #'
    #' dataset <- sequence_data$new("my_dataset")
    #' dataset$add_sequences(names)
    #' dataset$assign_sequence_taxonomy(names, taxonomies)
    #'
    assign_sequence_taxonomy = function(names, taxonomies) {
      self$data$assign_sequence_taxonomy(names, taxonomies)
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
    #' Remove all sequences from dataset
    clear = function() {
      self$data$clear()
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
    #' Get data frame containing sequence otu assignments
    #' @examples
    #'   dataset <- sequence_data$new("my_dataset")
    #'   seq_ids <- c("seq1", "seq2", "seq4", "seq3", "seq6", "seq5")
    #'   sequence_abundances <- c(10, 100, 1, 500, 25, 80)
    #'   otu_ids <- c("otu1", "otu1", "otu1", "otu2", "otu2", "otu3")
    #'   dataset$assign_otus(otu_ids, sequence_abundances, seq_ids = seq_ids)
    #'
    #'   # (list) otus would look like:
    #'   # otu1             otu2        otu3
    #'   # seq1,seq2,seq4   seq3,seq6   seq5
    #'
    #'   list <- dataset$get_list()
    #'
    #'   #  list$otu_id  list$seq_id
    #'   #  otu1          seq1
    #'   #  otu1          seq2
    #'   #  otu1          seq4
    #'   #  otu2          seq3
    #'   #  otu2          seq6
    #'   #  otu3          seq5
    #'
    #' @return data.frame
    get_list = function() {
      self$data$get_list()
    },

    #' @description
    #' Get the number of otus in the dataset
    #' @return An integer
    get_num_otus = function() {
      self$data$num_otus
    },

    #' @description
    #' Get number of samples in the dataset
    #' @return An integer
    get_num_samples = function() {
      self$data$num_samples
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
      self$data$num_treatments
    },

    #' @description
    #' Get data.frame containing the oligo data
    #' tag, oligo, diffs, oligo(optional), diffs(optional), sample(optional)
    #' @return data.frame
    get_oligos = function() {
      # TODO
    },

    #' @description
    #' Get data frame containing sequence otu assignments
    #' @examples
    #'   dataset <- sequence_data$new("my_dataset")
    #'   otu_ids <- c("otu1", "otu2", "otu3")
    #'   otu_abundances <- c(111, 525, 80)
    #'   dataset$assign_otus(otu_ids, otu_abundances)
    #'
    #'   # (rabund) otus would look like:
    #'   # otu1  otu2  otu3
    #'   # 110   525   80
    #'
    #'   rabund <- dataset$get_rabund()
    #'
    #'   #  rabund$otu_id  rabund$abundance
    #'   #  otu1           111
    #'   #  otu2           525
    #'   #  otu3           80
    #'
    #' @return data.frame
    get_rabund = function() {
      self$data$get_rabund()
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
      if (self$data$num_samples != 0) {
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
    #' Get report containing the scrapped sequences / otus -
    #' ids, trash_codes
    #' @return list of data.frames
    get_scrap_report = function() {
      results <- list()
      list_names <- c("sequence_scrap_report")
      scrap_sequence_report <- self$data$get_scrap_report("sequence")
      results[[1]] <- scrap_sequence_report
      if (self$data$num_otus != 0) {
        results[[2]] <- self$data$get_scrap_report("otu")
        list_names <- c(list_names, "otu_scrap_report")
      }
      names(results) <- list_names
      results
    },

    #' @description
    #' Get report containing the sequence report table -
    #' ids, starts, ends, lengths, ambigs, homopolymers, numns
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
    #' dataset$add_sequences(names)
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
        return()
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
      if ((length(results[[1]]$starts) == 8) && (!silent)) {
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
    #' Get data frame containing sequence otu assignments by sample
    #' @examples
    #'   dataset <- sequence_data$new("my_dataset")
    #'   otu_ids <- c("otu1", "otu1", "otu1", "otu2", "otu2", "otu3")
    #'   samples <- c("sample1", "sample2", "sample5",
    #'    "sample1", "sample3", "sample1")
    #'   sample_abundances <- c(10, 100, 1, 500, 25, 80)
    #'   dataset$assign_otus(otu_ids, sample_abundances, samples)
    #'
    #'   # (shared) otus would look like:
    #'   # label  sample   otu1   otu2   otu3
    #'   # 0.03   sample1  10     500    80
    #'   # 0.03   sample2  100    0      0
    #'   # 0.03   sample3  0      25     0
    #'   # 0.03   sample5  1      0      0
    #'
    #'   shared <- dataset$get_shared()
    #'
    #'   # shared$otu_id  shared$sample  shared$abundance
    #'   # otu1           sample1        10
    #'   # otu1           sample2        100
    #'   # otu1           sample5        1
    #'   # otu2           sample1        500
    #'   # otu2           sample3        25
    #'   # otu3           sample1        80
    #'
    #' @return data.frame
    get_shared = function() {
      self$data$get_shared()
    },

    #' @description
    #' Get data.frames containing classifications for sequences and OTUs.
    #' id, taxonomy, abundance(optional)
    #' These tables represent mothur's taxonomy and cons.taxonomy files
    #' @return List of data.frames
    get_taxonomy_tables = function() {
      # TODO
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
    }
  ),
  private = list(
    # Clear sequences from dataset
    finalize = function() {
      self$clear()
    }
  )
)
