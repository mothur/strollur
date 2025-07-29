#' @title read_mothur
#' @description
#' The read_mothur function will read various
#' \href{https://mothur.org/wiki/tags/#file_types}{file types} created by
#' mothur, and create a "sequence_data" object.
#'
#' To generate the various output files you can follow Pat's
#' \href{https://mothur.org/wiki/miseq_sop/}{Miseq example analysis}.
#' @param dataset_name A string containing a name for your dataset.
#' @param fasta filename, a FASTA formatted file containing sequence strings.
#' \href{https://mothur.org/wiki/fasta_file/}{fasta file}
#' @param count filename, a mothur
#' \href{https://mothur.org/wiki/count_file/}{count file}
#' @param design filename, a mothur
#' \href{https://mothur.org/wiki/design_file/}{design file}
#' @param taxonomy filename, a mothur
#' \href{https://mothur.org/wiki/taxonomy_file/}{taxonomy file}
#' @param list filename, a mothur
#' \href{https://mothur.org/wiki/list_file/}{list file}
#' @param shared filename, a mothur
#' \href{https://mothur.org/wiki/shared_file/}{shared file}
#' @param rabund filename, a mothur
#' \href{https://mothur.org/wiki/rabund_file/}{rabund file}
#' @param cons_taxonomy filename, a mothur consensus taxonomy file
#' \href{https://mothur.org/wiki/constaxonomy_file/}{constaxonomy file}
#' @param list_type string containing the type of list file. Options include:
#' 'otu', 'asv' or 'phylotype'. Default = 'otu'.
#' @param shared_type string containing the type of shared file. Options
#' include: 'otu', 'asv' or 'phylotype'. Default = 'otu'.
#' @param rabund_type string containing the type of rabund file. Options
#'  include: 'otu', 'asv' or 'phylotype'. Default = 'otu'.
#'
#' @note
#' \itemize{
#' \item \emph{consensus taxonomy}, The 'sequence_data' object will generate
#' consensus taxonomies for you based on the sequence taxonomy assignment. You
#' only need to provide the ".cons.taxonomy" file if you are not providing
#' sequence taxonomy assignments.
#' \item \emph{shared / rabund file}, The 'sequence_data' object will generate
#' shared and rabund data for you based on the otu assignment in the list file
#' and the count data. You only need to provide the ".shared" file if you are
#' not providing the list and count files.
#' }
#' @examples
#' # For dataset's including sequence data:
#'
#' dataset <- read_mothur(
#'   fasta = rdataset_example("final.fasta"),
#'   count = rdataset_example("final.count_table"),
#'   taxonomy = rdataset_example("final.taxonomy"),
#'   design = rdataset_example("mouse.time.design"),
#'   list = rdataset_example("final.opti_mcc.list"),
#'   dataset_name = "miseq_sop"
#' )
#'
#' # For dataset's with only otu data:
#'
#' dataset <- read_mothur(
#'   shared = rdataset_example("final.opti_mcc.shared"),
#'   cons_taxonomy = rdataset_example(
#'     "final.cons.taxonomy"
#'   ),
#'   design = rdataset_example("mouse.time.design"),
#'   dataset_name = "miseq_sop"
#' )
#'
#' @return A 'sequence_data' object
#' @export
read_mothur <- function(fasta = NULL, count = NULL,
                        taxonomy = NULL, list = NULL, design = NULL,
                        shared = NULL, rabund = NULL, cons_taxonomy = NULL,
                        dataset_name = "", list_type = "otu",
                        shared_type = "otu", rabund_type = "otu") {
  # create new blank dataset
  dataset <- sequence_data$new(name = dataset_name)

  # add sequence nucleotide strings
  if (!is.null(fasta)) {
    fasta_data <- read_fasta(fasta)
    dataset$add_sequences(fasta_data$names, fasta_data$sequences)
  }

  # add sequence abundance data
  if (!is.null(count)) {
    count_table <- read_mothur_count(count)

    # you did not add fasta seqs
    if (is.null(fasta)) {
      dataset$add_sequences(unique(count_table$id))
    }

    # if the count file include samples, add them
    if ("sample" %in% names(count_table)) {
      dataset$assign_sequence_abundance(
        count_table$id,
        count_table$abundance,
        count_table$sample
      )
    } else {
      dataset$assign_sequence_abundance(
        count_table$id,
        count_table$abundance
      )
    }
  }

  # add taxonomy data
  if (!is.null(taxonomy)) {
    df <- readr::read_table(
      file = taxonomy, col_names = FALSE,
      show_col_types = FALSE
    )
    dataset$assign_sequence_taxonomy(df[[1]], df[[2]])
  }

  # add sample / treatment assignments
  if (!is.null(design)) {
    df <- readr::read_table(
      file = design, col_names = TRUE,
      show_col_types = FALSE
    )
    dataset$assign_treatments(df[[1]], df[[2]])
  }

  # add sequence otu assignments
  if (!is.null(list) || !is.null(shared) || !is.null(rabund)) {
    # only one file should be used to assign the otus unless the types are
    # different. ie. list file for distance based otus, and file for
    # phylotype otus, and rabund for asv's.
    if (!is.null(list)) {
      df <- read_mothur_list(list)
      dataset$assign_bins(df$bin_id,
        seq_ids = df$seq_id,
        type = list_type
      )
    }

    if (!is.null(rabund)) {
      df <- read_mothur_rabund(rabund)
      dataset$assign_bins(df$bin_id, df$abundance, type = rabund_type)
    }

    if (!is.null(shared)) {
      df <- read_mothur_shared(shared)
      dataset$assign_bins(df$bin_id, df$abundance,
        df$sample,
        type = shared_type
      )
    }
  }

  if (!is.null(cons_taxonomy)) {
    df <- readr::read_table(
      file = cons_taxonomy, col_names = TRUE,
      show_col_types = FALSE
    )

    # remove size
    df <- df[, -c(2)]
    dataset$assign_bin_taxonomy(df[[1]], df[[2]])
  }

  dataset
}
