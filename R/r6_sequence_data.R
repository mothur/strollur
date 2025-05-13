#' @title sequence_dataset
#' @description 'sequence_data' is an R6 class that represents sequence
#'  FASTA and abundance data.
#'
#' @author Sarah Westcott, \email{swestcot@@umich.edu}
#'
#' @importFrom R6 R6Class
#' @importFrom methods new
#' @import cli
#' @import data.table
#' @export
sequence_data <- R6Class("sequence_data",
  public = list(

    #' @field dataset pointer to C++ class
    dataset = NULL,

    #' @description
    #' Create a new FASTA sequence dataset
    #' @param name Name of dataset
    #' @examples
    #'   dataset <- sequence_data$new("soil")
    #'
    #' @return A new `sequence_data` object.
    initialize = function(name) {
      self$dataset <- new(Dataset, name)
      invisible(self)
    },

    #' @description
    #' Get summary of sequences data
    print = function() {
      message(self$dataset$print())
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
        self$dataset$add_seqs(names, sequences, comments)

      invisible(self)
    },

    #' @description
    #' Remove all sequences from dataset
    clear = function() {
      self$dataset$clear()
      invisible(self)
    },

    #' @description
    #' Get List containing dataset
    #' @return List
    export = function() {
      return(self$dataset$export())
    },

    #' @description
    #' Get names of groups in the dataset
    #' @param name The name of the sequence you want groups for, optional
    #' @return A character vector
    get_groups = function(name = NULL) {
        if (is.null(name)) {
            name = ""
        }
      return(self$dataset$get_groups(name))
    },

    #' @description
    #' Get the number of sequences represented in samples
    #' @param group String, name of sample
    #' @return A vector of integers
    get_group_totals = function(group = NULL) {
      if (is.null(group)) {
        group = ""
      }
      return(self$dataset$get_group_totals(group))
    },

    #' @description
    #' Get names of sequences in the dataset
    #' @param group String, name of sample
    get_names = function(group = NULL) {
      if (is.null(group)) {
        group = ""
      }
      return(self$dataset$get_names(group))
    },

    #' @description
    #' Get number of groups in the dataset
    #' @return An integer
    get_num_groups = function() {
      return(self$dataset$num_groups)
    },

    #' @description
    #' Get the number of sequences in the dataset
    #' @param group The name of the group you want number of sequences for,
    #'  optional
    #' @return An integer
    get_num_seqs = function(group = NULL) {
      if (is.null(group)) {
        group = ""
      }
      return(self$dataset$get_total(group))
    },

    #' @description
    #' Get number of unique sequences in dataset
    #' @return An integer
    get_num_unique = function() {
      return(self$dataset$num_unique)
    },

    #' @description
    #' Get sample table returns data.table containing the 3 columns: name,
    #' group (optional) and abundance.
    #' @return data.table
    get_sample_table = function() {
      return(self$dataset$get_sequence_abundance_table())
    },


    #' @description
    #' Get vector containing FASTA nucleotide strings
    #' @param group String, name of sample
    get_seqs = function(group = NULL) {
      if (is.null(group)) {
        group = ""
      }

      return(self$dataset$get_seqs(group))
    },

    #' @description
    #' Determine if a group is present in the dataset
    #' @param group String, Name of sample
    #' @return Boolean
    has_group = function(group) {
      return(self$dataset$has_group(group))
    },

    #' @description
    #' Determine if the dataset is aligned and its alignment length
    #' @return bool
    is_aligned = function() {
       return (self$dataset$is_aligned)
    },

    #' @description
    #' Reinstates sequences removed for the reasons in trash_codes
    #' @param trash_codes vector containing reasons for sequences removal
    reinstate_seqs = function(trash_codes) {
      if (length(trash_codes) != 0) {
        self$dataset$reinstate_seqs(trash_codes)
      }

      invisible(self)
    },

    #' @description
    #' Remove sequences from dataset for cause
    #' @param names vector containing names of sequences to remove
    #' @param trash_codes vector containing reasons for sequences removal
    remove_seqs = function(names, trash_codes) {
      if (length(names) != length(trash_codes)) {
        super$abort_length_mismatch(
          "names", "trash_codes", length(names),
          length(trash_codes)
        )
      } else {
        self$dataset$remove_seqs(names, trash_codes)
      }
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
        self$dataset$assign_sample_abundance(names, groups, abundances)
      }

      invisible(self)
    },

    #' @description
    #' Set sequence data
    #' @param names a vector of sequence names
    #' @param sequences a vector of sequence data
    #' @param comments a vector of sequence comments, (optional)
    set_seqs = function(names, sequences, comments = NULL) {
        if (is.null(comments)) {
            comments <- rep("", length(names))
        }
        self$dataset$set_seqs(names, sequences, comments)

        invisible(self)
    }
  ),

  private = list(
    # Clear sequences from dataset
    finalize = function() {
      self$clear()
    }
  )
)
