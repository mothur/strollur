#' @title sequence_dataset
#' @description 'sequence_data' is an R6 class that represents sequence
#'  FASTA and abundance data.
#'
#' @author Sarah Westcott, \email{swestcot@@umich.edu}
#'
#' @importFrom R6 R6Class
#' @importFrom methods new
#' @importFrom parallelly, availableCores
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
        fasta_data <- read_fasta_file(fasta)
        self$add_seqs(fasta_data$names, fasta_data$sequences)
      }
      invisible(self)
    },

    #' @description
    #' Get summary of sequences data
    print = function() {
      message(self$data$print())
    },

    #' @description
    #' Add new sequence data
    #' @param names a vector of sequence names
    #' @param sequences a vector of sequence data
    #' @param comments a vector of sequence comments, (optional)
    #' @examples
    #'
    #'   dataset <- sequence_data$new("my_dataset")
    #'   fasta_data <- read_fasta_file(rdataset_example("test.fasta"))
    #'   dataset$add_seqs(fasta_data$names, fasta_data$sequences)
    #'
    add_seqs = function(names, sequences, comments = NULL) {
      if (is.null(comments)) {
        comments <- rep("", length(names))
      }
      self$data$add_seqs(names, sequences, comments)

      invisible(self)
    },

    #' @description
    #' Add group assignment data
    #' @param names a vector of sequence names
    #' @param groups a vector of group assignments
    #' @param abundances a vector of sample abundances
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
    #' groups <- c("sample2", "sample3", "sample4",
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
    #' dataset$assign_sample_abundance(names, groups, abundances)
    #'
    assign_sample_abundance = function(names = NULL, groups = NULL,
                                       abundances = NULL) {
      if (!is.null(names) && !is.null(groups)) {
        self$data$assign_sample_abundance(names, groups, abundances)
      }

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
    #' @return data.table
    get_accnos_summary = function() {
      # TODO
    },

    #' @description
    #' Get report containing align report table -
    #' search_scores, sim_scores, longest_insert
    #' @return data.table
    get_align_report = function() {
      # TODO
    },

    #' @description
    #' Get summary of the align report
    #' @return data.table
    get_align_summary = function() {
      # TODO
    },

    #' @description
    #' Get report containing contigs report table -
    #' ostarts, oends, olengths, mismatches, ee
    #' @return data.table
    get_contigs_report = function() {
      # TODO
    },

    #' @description
    #' Get summary of the contigs report
    #' @return data.table
    get_contigs_summary = function() {
      # TODO
    },

    #' @description
    #' Get count table returns data.table containing:
    #' ids, abundances, group(optional), treatment(optional)
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
    #' Get sparse column formatted distance table
    #' @return data.table
    get_dist = function() {
      # TODO
    },

    #' @description
    #' Get report containing the fasta report table -
    #' starts, ends, lengths, ambigs, homopolymers, numns
    #' @return data.table
    get_fasta_report = function() {
        fasta_report <- self$data$get_fasta_report()
        data.table(starts = fasta_report[[1]],
                   ends = fasta_report[[2]],
                   lengths = fasta_report[[3]],
                   ambigs = fasta_report[[4]],
                   homopolymers = fasta_report[[5]],
                   numns = fasta_report[[6]])
    },

    #' @description
    #' Get summary of the fasta report
    #' @param silent Default = FALSE, meaning print fasta summary
    #' @return data.table
    get_fasta_summary = function(silent = FALSE) {
        fasta_results <- self$data$get_fasta_summary()

        results_row_names <- c(
            "Minimum:", "2.5%-tile:", "25%-tile:",
            "Median:   ", "75%-tile:", "97.5%-tile:",
            "Maximum:", "Mean:      "
        )

        rownames(fasta_results) <- results_row_names

        # if you have summary results to print
        if ((length(fasta_results$starts) == 8) && (!silent)) {

            cat("\t\tstart\tend\tlength\tambigs\tpolymer\tnum_ns\tnumseqs\n")
            cat(
                paste(results_row_names[1], fasta_results$starts[1],
                      fasta_results$ends[1],
                      fasta_results$nbases[1], fasta_results$ambigs[1],
                      fasta_results$polymers[1], fasta_results$numns[1],
                      fasta_results$numseqs[1],
                      sep = "\t"
                ),
                "\n"
            )
            cat(
                paste(results_row_names[2], fasta_results$starts[2],
                      fasta_results$ends[2],
                      fasta_results$nbases[2], fasta_results$ambigs[2],
                      fasta_results$polymers[2], fasta_results$numns[2],
                      fasta_results$numseqs[2],
                      sep = "\t"
                ),
                "\n"
            )
            cat(
                paste(results_row_names[3], fasta_results$starts[3],
                      fasta_results$ends[3],
                      fasta_results$nbases[3], fasta_results$ambigs[3],
                      fasta_results$polymers[3], fasta_results$numns[3],
                      fasta_results$numseqs[3],
                      sep = "\t"
                ),
                "\n"
            )
            cat(
                paste(results_row_names[4], fasta_results$starts[4],
                      fasta_results$ends[4],
                      fasta_results$nbases[4], fasta_results$ambigs[4],
                      fasta_results$polymers[4], fasta_results$numns[4],
                      fasta_results$numseqs[4],
                      sep = "\t"
                ),
                "\n"
            )
            cat(
                paste(results_row_names[5], fasta_results$starts[5],
                      fasta_results$ends[5],
                      fasta_results$nbases[5], fasta_results$ambigs[5],
                      fasta_results$polymers[5], fasta_results$numns[5],
                      fasta_results$numseqs[5],
                      sep = "\t"
                ),
                "\n"
            )
            cat(
                paste(results_row_names[6], fasta_results$starts[6],
                      fasta_results$ends[6],
                      fasta_results$nbases[6], fasta_results$ambigs[6],
                      fasta_results$polymers[6], fasta_results$numns[6],
                      fasta_results$numseqs[6],
                      sep = "\t"
                ),
                "\n"
            )
            cat(
                paste(results_row_names[7], fasta_results$starts[7],
                      fasta_results$ends[7],
                      fasta_results$nbases[7], fasta_results$ambigs[7],
                      fasta_results$polymers[7], fasta_results$numns[7],
                      fasta_results$numseqs[7],
                      sep = "\t"
                ),
                "\n"
            )
            cat(paste(
                results_row_names[8], fasta_results$starts[8], fasta_results$ends[8],
                fasta_results$nbases[8], fasta_results$ambigs[8],
                fasta_results$polymers[8], fasta_results$numns[8],
                sep = "\t"
            ), "\n")
            cat("Unique seqs:\t", self$get_num_unique(), "\n")
            cat("Total seqs:\t", self$get_num_seqs(), "\n")
        }
        return (fasta_results)
    },

    #' @description
    #' Get names of groups in the dataset
    #' @param name The name of the sequence you want groups for, optional
    #' @return A character vector
    get_groups = function(name = NULL) {
      if (is.null(name)) {
        name <- ""
      }
      self$data$get_groups(name)
    },

    #' @description
    #' Get summary of count table. group, group_totals
    #' @return data.table
    get_group_summary = function() {
      if (is.null(group)) {
        group <- ""
      }
      self$data$get_group_totals(group)
    },

    #' @description
    #' Get names of sequences in the dataset
    #' @param group String, name of sample
    get_ids = function(group = NULL) {
      if (is.null(group)) {
        group <- ""
      }
      self$data$get_names(group)
    },

    #' @description
    #' Get number of groups in the dataset
    #' @return An integer
    get_num_groups = function() {
      self$data$num_groups
    },

    #' @description
    #' Get the number of sequences in the dataset
    #' @param group The name of the group you want number of sequences for,
    #'  optional
    #' @return An integer
    get_num_seqs = function(group = NULL) {
      if (is.null(group)) {
        group <- ""
      }
      self$data$get_total(group)
    },

    #' @description
    #' Get number of unique sequences in dataset
    #' @return An integer
    get_num_unique = function() {
      self$data$num_unique
    },

    #' @description
    #' Get data.table containing the oligo data
    #' tag, oligo, diffs, oligo(optional), diffs(optional), group(optional)
    #' @return data.table
    get_oligos = function() {
      # TODO
    },

    #' @description
    #' Get data.tables representing your OTU / ASV clusters.
    #' id, abundance, group(optional)
    #' These tables represent mothur's list, rabund, shared, relabund files
    #' @return List of data.tables
    get_otu_tables = function() {
      # TODO
    },

    #' @description
    #' Get vector containing FASTA nucleotide strings
    #' @param group String, name of sample
    get_seqs = function(group = NULL) {
      if (is.null(group)) {
        group <- ""
      }

      self$data$get_seqs(group)
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
    #' Determine if a group is present in the dataset
    #' @param group String, Name of sample
    #' @return Boolean
    has_group = function(group) {
      self$data$has_group(group)
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
