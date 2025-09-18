#' @title miseq_sop_example
#' @description
#' The miseq_sop_example function will create 'sequence_data' object using the
#' analysis files from the \href{https://mothur.org/wiki/miseq_sop/}{MiSeq_SOP}
#' example.
#'
#' @examples
#'
#' miseq_dataset <- miseq_sop_example()
#'
#' @return A 'sequence_data' object
#' @export
miseq_sop_example <- function() {
  dataset <- read_mothur(
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

  metadata <- readr::read_tsv(rdataset_example("mouse.dpw.metadata"),
    col_names = TRUE, show_col_types = FALSE
  )
  dataset$add_metadata(metadata)

  reference <- readr::read_csv(rdataset_example("references.csv"),
    col_names = TRUE, show_col_types = FALSE
  )
  dataset$add_references(reference = reference)

  dataset
}
