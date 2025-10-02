# tests import of dataset object

test_that("import - miseq_sop_example", {
  # full dataset
  miseq <- miseq_sop_example()

  expect_equal(miseq$get_num_bins("phylotype"), 63)
  expect_equal(miseq$get_num_sequences(TRUE), 2425)

  phylo_bins_to_remove <- c("Phylo01", "Phylo02")
  reasons_to_remove <- c("testing", "testing")

  # remove some bins to allow for filtering
  remove_bins(
    miseq$data, phylo_bins_to_remove,
    reasons_to_remove, "phylotype"
  )

  expect_equal(miseq$get_num_bins("phylotype"), 61)
  expect_equal(miseq$get_num_sequences(TRUE), 825)
  expect_equal(miseq$get_num_sequences(), 39177)

  exported_miseq <- miseq$export()

  expect_equal(sum(exported_miseq$sequence_data$include_sequence), 825)

  dataset_t <- import(exported_miseq)

  expect_equal(dataset_t$get_dataset_name(), miseq$get_dataset_name())
  expect_equal(dataset_t$get_num_sequences(TRUE), miseq$get_num_sequences(TRUE))
  expect_equal(dataset_t$get_num_sequences(), miseq$get_num_sequences())
  expect_equal(dataset_t$get_num_treatments(), miseq$get_num_treatments())
  expect_equal(dataset_t$get_num_samples(), miseq$get_num_samples())
  expect_equal(dataset_t$get_num_bins("otu"), miseq$get_num_bins("otu"))
  expect_equal(
    dataset_t$get_num_bins("phylotype"),
    miseq$get_num_bins("phylotype")
  )
  expect_equal(dataset_t$get_num_bins("asv"), miseq$get_num_bins("asv"))
  expect_equal(
    get_sample_totals(dataset_t$data),
    get_sample_totals(miseq$data)
  )
  expect_equal(
    get_treatment_totals(dataset_t$data),
    get_treatment_totals(miseq$data)
  )

  # only import bin data no sequences
  dataset_t <- import(exported_miseq, c("bin_data"))

  # abundances reflect the total abundance of the sequence in bin
  expect_equal(dataset_t$get_rabund("asv")$abundance[1:3], c(7436, 6285, 5207))
  # no list, since no sequences
  expect_equal(dataset_t$get_list("asv"), data.frame())
  expect_equal(dataset_t$get_sequence_taxonomy_report(), data.frame())
  # 303 bins x 6 tax levels
  expect_equal(nrow(dataset_t$get_bin_taxonomy_report()), 1818)

  data <- dataset$new()
  abunds <- c(1, 10, 100)
  bins <- c("otu1", "otu2", "otu3")
  data$add_alignment_report(
    readr::read_tsv(rdataset_example("alignment_data.tsv"),
      col_names = TRUE, show_col_types = FALSE
    ),
    "QueryName"
  )
  data$add_contigs_assembly_report(
    readr::read_tsv(rdataset_example("contigs_data.tsv"),
      col_names = TRUE, show_col_types = FALSE
    ), "Name"
  )
  data$assign_bins(bin_names = bins, abundances = abunds)

  expect_equal(data$get_num_bins(), 3)
  expect_equal(data$get_num_sequences(TRUE), 5)
  expect_equal(data$get_num_sequences(), 111)
  expect_equal(data$get_num_samples(), 0)
  expect_equal(data$get_num_treatments(), 0)

  table <- export(data)
  dataset2 <- import(table)

  expect_equal(dataset2$get_num_bins(), 3)
  expect_equal(dataset2$get_num_sequences(TRUE), 5)
  expect_equal(dataset2$get_num_sequences(), 111)
  expect_equal(dataset2$get_num_samples(), 0)
  expect_equal(dataset2$get_num_treatments(), 0)
})

test_that("import - no sequence data", {
  # just shared and constax dataset
  just_bins <- read_mothur(
    otu_shared = rdataset_example("final.opti_mcc.shared"),
    cons_taxonomy = rdataset_example("final.cons.taxonomy"),
    design = rdataset_example("mouse.time.design"),
    sample_tree = rdataset_example("final.opti_mcc.jclass.ave.tre"),
    dataset_name = "just_bins"
  )

  expect_equal(just_bins$get_num_bins("otu"), 531)
  expect_equal(just_bins$get_num_sequences(), 113963)

  table <- just_bins$export()

  dataset <- import(table)

  expect_equal(dataset$get_num_bins("otu"), 531)
  expect_equal(dataset$get_num_sequences(), 113963)

  expect_equal(dataset$get_dataset_name(), just_bins$get_dataset_name())

  expect_equal(dataset$get_num_sequences(), just_bins$get_num_sequences())
  expect_equal(dataset$get_num_treatments(), just_bins$get_num_treatments())
  expect_equal(dataset$get_num_samples(), just_bins$get_num_samples())
  expect_equal(length(dataset$get_sequence_names()), 0)

  expect_error(import(table, c("sequence_data")))
  expect_error(import(table, c("metadata")))
  expect_error(import(table, c("references")))
  expect_error(import(table, c("alignment_report")))
  expect_error(import(table, c("contigs_assembly_report")))
  expect_error(import(table, c("sequence_tree")))

  attr(table, "rdataset_version") <- "0.2.1"
  expect_error(import(table))
})

test_that("import - no bin data", {
  # just shared and constax dataset
  just_seqs <- read_mothur(
    fasta = rdataset_example("final.fasta"),
    count = rdataset_example("final.count_table"),
    taxonomy = rdataset_example("final.taxonomy"),
    design = rdataset_example("mouse.time.design"),
    sequence_tree = rdataset_example("final.phylip.tre"),
    dataset_name = "just_seqs"
  )

  expect_equal(just_seqs$get_num_bins("otu"), 0)
  expect_equal(just_seqs$get_num_sequences(), 113963)
  expect_equal(just_seqs$get_num_sequences(TRUE), 2425)

  table <- just_seqs$export()

  dataset <- import(table)

  expect_equal(dataset$get_num_bins("otu"), 0)
  expect_equal(dataset$get_num_sequences(), 113963)
  expect_equal(just_seqs$get_num_sequences(TRUE), 2425)

  expect_equal(dataset$get_dataset_name(), just_seqs$get_dataset_name())

  expect_equal(dataset$get_num_sequences(), just_seqs$get_num_sequences())
  expect_equal(dataset$get_num_treatments(), just_seqs$get_num_treatments())
  expect_equal(dataset$get_num_samples(), just_seqs$get_num_samples())
  expect_equal(length(dataset$get_sequence_names()), 2425)

  expect_error(import(table, c("bin_data")))
  expect_error(import(table, c("sample_tree")))
})

test_that("import - errors and warnings", {
  just_bins <- read_mothur(
    otu_shared = rdataset_example("final.opti_mcc.shared"),
    cons_taxonomy = rdataset_example("final.cons.taxonomy"),
    design = rdataset_example("mouse.time.design"),
    sample_tree = rdataset_example("final.opti_mcc.jclass.ave.tre"),
    dataset_name = "just_bins"
  )

  expect_equal(just_bins$get_num_bins("otu"), 531)
  expect_equal(just_bins$get_num_sequences(), 113963)

  table <- just_bins$export()

  dataset <- import(table)

  expect_equal(dataset$get_num_bins("otu"), 531)
  expect_equal(dataset$get_num_sequences(), 113963)

  expect_equal(dataset$get_dataset_name(), just_bins$get_dataset_name())

  expect_equal(dataset$get_num_sequences(), just_bins$get_num_sequences())
  expect_equal(dataset$get_num_treatments(), just_bins$get_num_treatments())
  expect_equal(dataset$get_num_samples(), just_bins$get_num_samples())
  expect_equal(length(dataset$get_sequence_names()), 0)
})
