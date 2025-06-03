#' @title sequence_dataset
#' @description 'sequence_data' is an R6 class that represents sequence
#'  FASTA and abundance data.
#'
#' @author Sarah Westcott, \email{swestcot@@umich.edu}
#'
#' @importFrom R6 R6Class
#' @importFrom methods new
#' @importFrom parallelly, availableCores
#' @importFrom waldo compare
#' @import cli
#' @import data.table
#' @export
sequence_data <- R6Class("sequence_data",
  public = list(

    #' @field data pointer to C++ class
    data = NULL,

    #' @description
    #' Create a new FASTA sequence dataset
    #' @param name Name of dataset
    #' @param fasta, String name of FASTA file
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
        self$add_seqs(fasta_data$names, fasta_data$sequences)
      }
      invisible(self)
    },

    #' @description
    #' Get summary of sequence data
    print = function() {
      self$get_accnos_summary()
      cat("\n\n")
      self$get_sequence_summary()
      cat("\n\n")
      self$get_sample_summary()
      cat(
        paste("\nNumber of unique seqs:", self$get_num_sequences(TRUE)),
        "\n"
      )
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
    #'   sequence_data <- read_fasta(rdataset_example("test.fasta"))
    #'   dataset$add_seqs(sequence_data$names, sequence_data$sequences)
    #'
    add_seqs = function(names, sequences, comments = NULL) {
      if (is.null(comments)) {
        comments <- rep("", length(names))
      }
      self$data$add_seqs(names, sequences, comments)

      invisible(self)
    },

    #' @description
    #' Add sequence abundance data with optional sample and treatment assignments
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
    #' names <- c("seq1", "seq1", "seq1",
    #'           "seq2", "seq2", "seq2",
    #'           "seq3", "seq3",
    #'           "seq4")
    #' samples <- c("sample2", "sample3", "sample4",
    #'            "sample2", "sample3", "sample4",
    #'            "sample2", "sample3",
    #'            "sample4")
    #' abundances <- c(250, 400, 500,
    #'                25, 40, 50,
    #'                25, 25,
    #'                4)
    #'
    #' dataset <- sequence_data$new("mydata")
    #' unique_names <- unique(names)
    #' sequences <- c("ATGGGCT", "..TG--ACCGT..", "..GGuatgc..", "..GGTAC-T..")
    #' dataset$add_seqs(unique_names, sequences)
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
    #' Get report containing eliminated sequence names and trash_codes
    #' @return data.table
    get_accnos_report = function() {
      # TODO
    },

    #' @description
    #' Get summary of eliminated sequences
    #' @param silent Default = FALSE, meaning print accnos summary
    #' @return data.table
    get_accnos_summary = function(silent = FALSE) {
      # TODO
    },

    #' @description
    #' Get report containing align report table -
    #' ids, search_scores, sim_scores, longest_insert
    #' @return data.table
    get_align_report = function() {
      if (self$data$has_align_data) {
        return(self$data$get_align_report())
      }
      data.table()
    },

    #' @description
    #' Get report containing contigs report table -
    #' ids, lengths, ostarts, oends, olengths, mismatches, numns, ee
    #' @return data.table
    get_contigs_report = function() {
      if (self$data$has_contigs_data) {
        return(self$data$get_contigs_report())
      }
      data.table()
    },

    #' @description
    #' Get count table returns data.table containing:
    #' ids, abundances, sample(optional), treatment(optional)
    #' This table represents mothur's count and design files.
    #' @return data.table
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
    #' Get report containing the sequence report table -
    #' ids, starts, ends, lengths, ambigs, homopolymers, numns
    #' @return data.table
    get_sequence_report = function() {
      sequence_report <- self$data$get_sequence_report()
    },

    #' @description
    #' Get summary of the sequence reports
    #' @param silent Default = FALSE, meaning print summaries
    #' @return list of data.tables
    get_sequence_summary = function(silent = FALSE) {
      results <- self$data$get_sequence_summary()

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

      return(results)
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
              cat(paste(treatment_names[i], treatment_totals[i], sep = "\t"), "\n")
            }
          }
        }

        if (self$get_num_treatments() != 0) {
          treatment_names <- self$get_treatments()
          treatment_totals <- self$data$get_treatment_totals()
          return(list(
            data.table(sample = sample_names, total = sample_totals),
            data.table(treatment = treatment_names, total = treatment_totals)
          ))
        } else {
          return(list(data.table(sample = sample_names, total = sample_totals)))
        }
      } else {
        cli::cli_alert("Your dataset does not include sample data, ignoring.")
      }
      list()
    },

    #' @description
    #' Get names of sequences in the dataset
    #' @param sample String, name of sample
    get_ids = function(sample = NULL) {
      if (is.null(sample)) {
        sample <- ""
      }
      self$data$get_names(sample)
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
    #' Get data.table containing the oligo data
    #' tag, oligo, diffs, oligo(optional), diffs(optional), sample(optional)
    #' @return data.table
    get_oligos = function() {
      # TODO
    },

    #' @description
    #' Get data.tables representing your OTU / ASV clusters.
    #' id, abundance, sample(optional)
    #' These tables represent mothur's list, rabund, shared, relabund files
    #' @return List of data.tables
    get_otu_tables = function() {
      # TODO
    },

    #' @description
    #' Get vector containing FASTA nucleotide strings
    #' @param sample String, name of sample
    get_seqs = function(sample = NULL) {
      if (is.null(sample)) {
        sample <- ""
      }

      self$data$get_seqs(sample)
    },

    #' @description
    #' Get data.tables containing classifications for sequences and OTUs.
    #' id, taxonomy, abundance(optional)
    #' These tables represent mothur's taxonomy and cons.taxonomy files
    #' @return List of data.tables
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
