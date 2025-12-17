#' @title miseq_sop_example
#' @description
#' The miseq_sop_example function will create 'dataset' object using the
#' analysis files from the \href{https://mothur.org/wiki/miseq_sop/}{MiSeq_SOP}
#' example.
#'
#' @examples
#'
#' miseq <- miseq_sop_example()
#'
#' @return A 'dataset' object
#' @export
miseq_sop_example <- function() {
  data <- read_mothur(
    fasta = rdataset_example("final.fasta"),
    count = rdataset_example("final.count_table"),
    taxonomy = rdataset_example("final.taxonomy"),
    cons_taxonomy = rdataset_example("final.cons.taxonomy"),
    design = rdataset_example("mouse.time.design"),
    otu_list = rdataset_example("final.opti_mcc.list"),
    asv_list = rdataset_example("final.asv.list"),
    phylo_list = rdataset_example("final.tx.list"),
    sample_tree = rdataset_example("final.opti_mcc.jclass.ave.tre"),
    sequence_tree = rdataset_example("final.phylip.tre"),
    dataset_name = "miseq_sop"
  )

  representative_seqs <- readr::read_tsv(
    rdataset_example("otu_representative_sequences.tsv"),
    col_names = TRUE, show_col_types = FALSE
  )
  assign(
    data = data, table = representative_seqs,
    type = "bin_representatives"
  )

  metadata <- readr::read_tsv(rdataset_example("mouse.dpw.metadata"),
    col_names = TRUE, show_col_types = FALSE
  )
  xdev_add_report(data, metadata, "metadata")

  reference <- readr::read_csv(rdataset_example("references.csv"),
    col_names = TRUE, show_col_types = FALSE
  )
  xdev_add_references(data, reference)

  contigs_report <- readr::read_tsv(rdataset_example("final.contigs_report"),
    col_names = TRUE, show_col_types = FALSE
  )
  add(
    data = data, table = contigs_report, type = "reports",
    report_type = "contigs_report",
    table_names = list(sequence_name = "Name")
  )

  data
}
