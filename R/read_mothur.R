#' @title read_mothur
#' @description
#' The read_mothur function will read various
#' \href{https://mothur.org/wiki/tags/#file_types}{file types} created by
#' mothur, and create a 'dataset' object.
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
#' \href{https://mothur.org/wiki/taxonomy_file/}{taxonomy file}, created by
#' \href{https://mothur.org/wiki/classify.seqs/}{classify.seqs}
#' @param otu_list filename, a mothur
#' \href{https://mothur.org/wiki/list_file/}{list file} containing otu bin
#'  assignments. The otu_list file is created by
#' \href{https://mothur.org/wiki/cluster/}{cluster},
#' \href{https://mothur.org/wiki/cluster.split/}{cluster.split}, and
#' \href{https://mothur.org/wiki/cluster.fit/}{cluster.fit}
#' @param asv_list filename, a mothur
#' \href{https://mothur.org/wiki/list_file/}{list file} containing asv bin
#'  assignments. The asv_list file is created by
#' \href{https://mothur.org/wiki/cluster/}{cluster} using the 'unique' method.
#' @param phylo_list filename, a mothur
#' \href{https://mothur.org/wiki/list_file/}{list file} containing phylotype bin
#'  assignments. The phylo_list file is created by
#' \href{https://mothur.org/wiki/phylotype/}{phylotype}.
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
#' \href{https://mothur.org/wiki/constaxonomy_file/}{constaxonomy file}. The
#' cons_taxonomy file is created by
#' \href{https://mothur.org/wiki/classify.otu/}{classify.otu}.
#' @param sample_tree filename, a tree that relates samples.
#' The sample tree is created by
#' \href{https://mothur.org/wiki/tree.shared/}{tree.shared}. We recommend
#'  running tree.shared with subsample = true, and using the 'ave.tre' output
#'  for best results.
#' @param sequence_tree filename, a tree that relates sequences.
#' The sequence tree is created by
#' \href{https://mothur.org/wiki/clearcut/}{clearcut}. We DO NOT recommend
#'  using sequence trees. With the ever growing size of modern datasets,
#'  sequence tree can be difficult / impossible to build without hitting a
#'  memory limitation.
#' @note
#' \itemize{
#' \item \emph{consensus taxonomy}, The 'dataset' object will generate
#' consensus taxonomies for you based on the sequence taxonomy assignment. You
#' only need to provide the ".cons.taxonomy" file if you are not providing
#' sequence taxonomy assignments.
#' \item \emph{shared / rabund file}, The 'dataset' object will generate
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
#'   sample_tree = rdataset_example("final.opti_mcc.jclass.ave.tre"),
#'   dataset_name = "miseq_sop"
#' )
#'
#' # For dataset's with only otu data:
#'
#' data <- read_mothur(
#'   otu_shared = rdataset_example("final.opti_mcc.shared"),
#'   cons_taxonomy = rdataset_example(
#'     "final.cons.taxonomy"
#'   ),
#'   design = rdataset_example("mouse.time.design"),
#'   sample_tree = rdataset_example("final.opti_mcc.jclass.ave.tre"),
#'   dataset_name = "miseq_sop"
#' )
#'
#' @return A 'dataset' object
#' @export
read_mothur <- function(fasta = NULL, count = NULL,
                        taxonomy = NULL, otu_list = NULL, asv_list = NULL,
                        phylo_list = NULL, design = NULL, cons_taxonomy = NULL,
                        otu_shared = NULL, asv_shared = NULL,
                        phylo_shared = NULL, sample_tree = NULL,
                        sequence_tree = NULL, dataset_name = "") {
  # create new blank dataset
  data <- dataset$new(name = dataset_name)

  # add sequence nucleotide strings
  if (!is.null(fasta)) {
    fasta_data <- read_fasta(fasta)
    data$add_sequences(fasta_data)
  }

  # add sequence abundance data
  if (!is.null(count)) {
    count_table <- read_mothur_count(count)

    # you did not add fasta seqs
    if (is.null(fasta)) {
      data$add_sequences(sequence_names = unique(count_table$id))
    }

    # if the count file include samples, add them
    if ("sample" %in% names(count_table)) {
      data$assign_sequence_abundance(
        data = NULL,
        count_table$id,
        count_table$abundance,
        count_table$sample
      )
    } else {
      data$assign_sequence_abundance(
        data = NULL,
        count_table$id,
        count_table$abundance
      )
    }
  }

  # add taxonomy data
  if (!is.null(taxonomy)) {
    df <- read_mothur_taxonomy(taxonomy)
    data$assign_sequence_taxonomy(df)
  }

  # add sequence otu assignments
  if (!is.null(otu_list)) {
    df <- read_mothur_list(otu_list)
    data$assign_bins(df, type = "otu")
  }

  if (!is.null(otu_shared)) {
    df <- read_mothur_shared(otu_shared)
    data$assign_bins(df, type = "otu")
  }

  if (!is.null(asv_list)) {
    df <- read_mothur_list(asv_list)
    data$assign_bins(df, type = "asv")
  }

  if (!is.null(asv_shared)) {
    df <- read_mothur_shared(asv_shared)
    data$assign_bins(df, type = "asv")
  }

  if (!is.null(phylo_list)) {
    df <- read_mothur_list(phylo_list)
    data$assign_bins(df, type = "phylotype")
  }

  if (!is.null(phylo_shared)) {
    df <- read_mothur_shared(phylo_shared)
    data$assign_bins(df, type = "phylotype")
  }

  # add sample / treatment assignments
  if (!is.null(design)) {
    df <- readr::read_table(
      file = design, col_names = TRUE,
      show_col_types = FALSE
    )
    data$assign_treatments(data = NULL, df[[1]], df[[2]])
  }

  if (!is.null(cons_taxonomy)) {
    df <- read_mothur_cons_taxonomy(cons_taxonomy)
    data$assign_bin_taxonomy(df)
  }

  if (!is.null(sample_tree)) {
    tree <- ape::read.tree(sample_tree)
    data$add_sample_tree(tree)
  }

  if (!is.null(sequence_tree)) {
    tree <- ape::read.tree(sequence_tree)
    data$add_sequence_tree(tree)
  }

  data
}
