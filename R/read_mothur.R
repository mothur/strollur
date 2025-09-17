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
#' @param otu_list filename, a mothur
#' \href{https://mothur.org/wiki/list_file/}{list file} containing otu bin
#'  assignments.
#' @param asv_list filename, a mothur
#' \href{https://mothur.org/wiki/list_file/}{list file} containing asv bin
#'  assignments.
#' @param phylo_list filename, a mothur
#' \href{https://mothur.org/wiki/list_file/}{list file} containing phylotype bin
#'  assignments.
#' @param otu_shared filename, a mothur
#' \href{https://mothur.org/wiki/shared_file/}{shared file} containing otu bin
#' sample abundance assignments.
#' @param asv_shared filename, a mothur
#' \href{https://mothur.org/wiki/shared_file/}{shared file} containing asv bin
#' sample abundance assignments.
#' @param phylo_shared filename, a mothur
#' \href{https://mothur.org/wiki/shared_file/}{shared file} containing phylotype
#'  bin sample abundance assignments.
#' @param cons_taxonomy filename, a mothur consensus taxonomy file
#' \href{https://mothur.org/wiki/constaxonomy_file/}{constaxonomy file}

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
#'   otu_list = rdataset_example("final.opti_mcc.list"),
#'   asv_list = rdataset_example("final.asv.list"),
#'   phylo_list = rdataset_example("final.tx.list"),
#'   dataset_name = "miseq_sop"
#' )
#'
#' # For dataset's with only otu data:
#'
#' dataset <- read_mothur(
#'   otu_shared = rdataset_example("final.opti_mcc.shared"),
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
                        taxonomy = NULL, otu_list = NULL, asv_list = NULL,
                        phylo_list = NULL, design = NULL, cons_taxonomy = NULL,
                        otu_shared = NULL, asv_shared = NULL,
                        phylo_shared = NULL, dataset_name = "") {
  # create new blank dataset
  dataset <- sequence_data$new(name = dataset_name)

  # add sequence nucleotide strings
  if (!is.null(fasta)) {
    fasta_data <- read_fasta(fasta)
    dataset$add_sequences(fasta_data)
  }

  # add sequence abundance data
  if (!is.null(count)) {
    count_table <- read_mothur_count(count)

    # you did not add fasta seqs
    if (is.null(fasta)) {
      dataset$add_sequences(sequence_names = unique(count_table$id))
    }

    # if the count file include samples, add them
    if ("sample" %in% names(count_table)) {
      dataset$assign_sequence_abundance(
        data = NULL,
        count_table$id,
        count_table$abundance,
        count_table$sample
      )
    } else {
      dataset$assign_sequence_abundance(
        data = NULL,
        count_table$id,
        count_table$abundance
      )
    }
  }

  # add taxonomy data
  if (!is.null(taxonomy)) {
    df <- read_mothur_taxonomy(taxonomy)
    dataset$assign_sequence_taxonomy(df)
  }

  # add sequence otu assignments
  if (!is.null(otu_list)) {
    df <- read_mothur_list(otu_list)
    dataset$assign_bins(df, type = "otu")
  }

  if (!is.null(otu_shared)) {
    df <- read_mothur_shared(otu_shared)
    dataset$assign_bins(df, type = "otu")
  }

  if (!is.null(asv_list)) {
    df <- read_mothur_list(asv_list)
    dataset$assign_bins(df, type = "asv")
  }

  if (!is.null(asv_shared)) {
    df <- read_mothur_shared(asv_shared)
    dataset$assign_bins(df, type = "asv")
  }

  if (!is.null(phylo_list)) {
    df <- read_mothur_list(phylo_list)
    dataset$assign_bins(df, type = "phylotype")
  }

  if (!is.null(phylo_shared)) {
    df <- read_mothur_shared(phylo_shared)
    dataset$assign_bins(df, type = "phylotype")
  }

  # add sample / treatment assignments
  if (!is.null(design)) {
    df <- readr::read_table(
      file = design, col_names = TRUE,
      show_col_types = FALSE
    )
    dataset$assign_treatments(data = NULL, df[[1]], df[[2]])
  }

  if (!is.null(cons_taxonomy)) {
    df <- read_mothur_cons_taxonomy(cons_taxonomy)
    dataset$assign_bin_taxonomy(df)
  }

  dataset
}
