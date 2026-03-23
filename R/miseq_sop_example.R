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

  representative_seqs <- read_tsv(
    strollur_example("otu_representative_sequences.tsv"),
    col_names = TRUE, show_col_types = FALSE
  )
  assign(
    data = data, table = representative_seqs,
    type = "bin_representatives"
  )

  metadata <- read_tsv(strollur_example("mouse.dpw.metadata"),
    col_names = TRUE, show_col_types = FALSE
  )
  xdev_add_report(data, metadata, "metadata")

  reference <- read_csv(strollur_example("references.csv"),
    col_names = TRUE, show_col_types = FALSE
  )
  xdev_add_references(data, reference)

  contigs_report <- read_tsv(strollur_example("final.contigs_report.gz"),
    col_names = TRUE, show_col_types = FALSE
  )
  add(
    data = data, table = contigs_report, type = "reports",
    report_type = "contigs_report",
    table_names = list(sequence_name = "Name")
  )

  data
}
