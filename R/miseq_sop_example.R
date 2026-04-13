#' @title Example strollur object
#' @description
#' The miseq_sop_example function will create 'strollur' object using the
#' analysis files from the \href{https://mothur.org/wiki/miseq_sop/}{MiSeq_SOP}
#' example.
#'
#' @examples
#'
#' miseq <- miseq_sop_example()
#'
#' @return A 'strollur' object
#' @export
miseq_sop_example <- function() {
  data <- read_mothur(
    fasta = strollur_example("final.fasta.gz"),
    count = strollur_example("final.count_table.gz"),
    taxonomy = strollur_example("final.taxonomy.gz"),
    cons_taxonomy = strollur_example("final.cons.taxonomy"),
    design = strollur_example("mouse.time.design"),
    otu_list = strollur_example("final.opti_mcc.list.gz"),
    asv_list = strollur_example("final.asv.list.gz"),
    phylo_list = strollur_example("final.tx.list.gz"),
    sample_tree = strollur_example("final.opti_mcc.jclass.ave.tre"),
    sequence_tree = strollur_example("final.phylip.tre.gz"),
    dataset_name = "miseq_sop"
  )

  representative_seqs <- readRDS(strollur_example(
    "miseq_representative_sequences.rds"
  ))

  xdev_assign_bin_representative_sequences(
    data = data, table = representative_seqs,
    bin_type = "otu"
  )

  metadata <- readRDS(strollur_example("miseq_metadata.rds"))
  xdev_add_report(data = data, table = metadata, type = "metadata")

  reference <- read_csv(strollur_example("references.csv"),
    col_names = TRUE, show_col_types = FALSE
  )
  xdev_add_references(data = data, table = reference)

  contigs_report <- readRDS(strollur_example("miseq_contigs_report.rds"))
  xdev_add_report(
    data = data, table = contigs_report, type = "contigs_report",
    sequence_name = "Name"
  )

  data
}
