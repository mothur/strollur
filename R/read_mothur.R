#' @title Create a
#'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#'   from mothur outputs
#' @name read_mothur
#' @rdname read_mothur
#' @description
#' The read_mothur function reads various
#' \href{https://mothur.org/wiki/tags/#file_types}{file types} created by
#' mothur, and creates a `strollur` object.
#'
#' To generate the various input files you can follow Pat's
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
#' \item \emph{consensus taxonomy}, The `strollur` object will generate
#' consensus taxonomies for you based on the sequence taxonomy assignment. You
#' only need to provide the ".cons.taxonomy" file if you are not providing
#' sequence taxonomy assignments.
#' \item \emph{shared / rabund file}, The `strollur` object will generate
#' shared and rabund data for you based on the otu assignment in the list file
#' and the count data. You only need to provide the ".shared" file if you are
#' not providing the list and count files.
#' }
#' @examples
#' # For dataset's including sequence data:
#'
#' data <- read_mothur(
#'   fasta = strollur_example("final.fasta.gz"),
#'   count = strollur_example("final.count_table.gz"),
#'   taxonomy = strollur_example("final.taxonomy.gz"),
#'   design = strollur_example("mouse.time.design"),
#'   otu_list = strollur_example("final.opti_mcc.list.gz"),
#'   asv_list = strollur_example("final.asv.list.gz"),
#'   phylo_list = strollur_example("final.tx.list.gz"),
#'   sample_tree = strollur_example("final.opti_mcc.jclass.ave.tre"),
#'   dataset_name = "miseq_sop"
#' )
#'
#' # For dataset's with only otu data:
#'
#' data <- read_mothur(
#'   otu_shared = strollur_example("final.opti_mcc.shared"),
#'   cons_taxonomy = strollur_example(
#'     "final.cons.taxonomy"
#'   ),
#'   design = strollur_example("mouse.time.design"),
#'   sample_tree = strollur_example("final.opti_mcc.jclass.ave.tre"),
#'   dataset_name = "miseq_sop"
#' )
#'
#' @return A
#'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#' @export
read_mothur <- function(fasta = NULL, count = NULL,
                        taxonomy = NULL, otu_list = NULL, asv_list = NULL,
                        phylo_list = NULL, design = NULL, cons_taxonomy = NULL,
                        otu_shared = NULL, asv_shared = NULL,
                        phylo_shared = NULL, sample_tree = NULL,
                        sequence_tree = NULL, dataset_name = "") {
  if (!is.null(otu_list) && !is.null(otu_shared)) {
    cli_abort("You can provide a list or shared file, not both.")
  }

  if (!is.null(asv_list) && !is.null(asv_shared)) {
    cli_abort("You can provide a list or shared file, not both.")
  }

  if (!is.null(phylo_list) && !is.null(phylo_shared)) {
    cli_abort("You can provide a list or shared file, not both.")
  }

  # create new blank dataset
  data <- strollur$new(dataset_name)

  # add sequence nucleotide strings
  if (!is.null(fasta)) {
    fasta_data <- read_fasta(fasta)
    xdev_add_sequences(data, fasta_data)
  }

  # add sequence abundance data
  if (!is.null(count)) {
    count_table <- read_mothur_count(count)

    # you did not add fasta seqs
    if (is.null(fasta)) {
      xdev_add_sequences(data, data.frame(
        sequence_names =
          unique(count_table$sequence_names)
      ))
    }

    xdev_assign_sequence_abundance(
      data,
      count_table
    )
  }

  # add taxonomy data
  if (!is.null(taxonomy)) {
    df <- read_mothur_taxonomy(taxonomy)
    xdev_assign_sequence_taxonomy(data, df)
  }

  # add sequence otu assignments
  if (!is.null(otu_list)) {
    df <- read_mothur_list(otu_list)
    assign(data = data, table = df, type = "bins", bin_type = "otu")
  }

  if (!is.null(otu_shared)) {
    df <- read_mothur_shared(otu_shared)
    assign(data = data, table = df, type = "bins", bin_type = "otu")
  }

  if (!is.null(asv_list)) {
    df <- read_mothur_list(asv_list)
    assign(data = data, table = df, type = "bins", bin_type = "asv")
  }

  if (!is.null(asv_shared)) {
    df <- read_mothur_shared(asv_shared)
    assign(data = data, table = df, type = "bins", bin_type = "asv")
  }

  if (!is.null(phylo_list)) {
    df <- read_mothur_list(phylo_list)
    assign(data = data, table = df, type = "bins", bin_type = "phylotype")
  }

  if (!is.null(phylo_shared)) {
    df <- read_mothur_shared(phylo_shared)
    assign(data = data, table = df, type = "bins", bin_type = "phylotype")
  }

  # add sample / treatment assignments
  if (!is.null(design)) {
    df <- readr::read_table(
      file = design, col_names = TRUE,
      show_col_types = FALSE
    )
    assign(data = data, table = df, type = "treatments")
  }

  if (!is.null(cons_taxonomy)) {
    df <- read_mothur_cons_taxonomy(cons_taxonomy)
    assign(data = data, table = df, type = "bin_taxonomy", bin_type = "otu")
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
