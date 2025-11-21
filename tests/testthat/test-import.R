# tests import of dataset object

test_that("import - miseq_sop_example", {
  # full dataset
  miseq <- miseq_sop_example()

  expect_equal(miseq$get_num_bins("phylotype"), 63)
  expect_equal(miseq$get_num_sequences(TRUE), 2425)

  phylo_bins_to_remove <- c("Phylo01", "Phylo02")
  reasons_to_remove <- c("testing", "testing")

  # remove some bins to allow for filtering
  xdev_remove_bins(
    miseq, phylo_bins_to_remove,
    reasons_to_remove, "phylotype"
  )

  expect_equal(miseq$get_num_bins("phylotype"), 61)
  expect_equal(miseq$get_num_sequences(TRUE), 825)
  expect_equal(miseq$get_num_sequences(), 39177)

  exported_miseq <- export_dataset(miseq)

  expect_equal(sum(exported_miseq$sequence_data$include_sequence), 825)

  dataset_t <- import_dataset(exported_miseq)

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
    get_sample_totals(dataset_t),
    get_sample_totals(miseq)
  )
  expect_equal(
    get_treatment_totals(dataset_t),
    get_treatment_totals(miseq)
  )

  dfd <- get_bin_representative_sequences(dataset_t)
  dfm <- get_bin_representative_sequences(miseq)

  for (otu in dfm[[1]]) {
    expect_equal(
      dfd %>% filter(otu_names == otu),
      dfm %>% filter(otu_names == otu)
    )
  }

  # only import bin data no sequences
  dataset_t <- import_dataset(exported_miseq, c("bin_data"))

  # abundances reflect the total abundance of the sequence in bin
  expect_equal(get_rabund(dataset_t, "asv")$abundance[1:3], c(7436, 6285, 5207))
  # no list, since no sequences
  expect_equal(get_list(dataset_t, "asv"), data.frame())
  expect_equal(report(dataset_t, "sequence_taxonomy"), data.frame())
  # 303 bins x 6 tax levels
  expect_equal(nrow(report(dataset_t, "bin_taxonomy")), 1818)

  data <- dataset$new()
  abunds <- c(1, 10, 100)
  bins <- c("otu1", "otu2", "otu3")
  add_report(
    data,
    readr::read_tsv(rdataset_example("alignment_data.tsv"),
      col_names = TRUE, show_col_types = FALSE
    ), "alignment_report",
    "QueryName"
  )
  add_report(
    data,
    readr::read_tsv(rdataset_example("contigs_data.tsv"),
      col_names = TRUE, show_col_types = FALSE
    ), "contigs_report", "Name"
  )
  assign_bins(data, data.frame(bin_names = bins, abundances = abunds))

  expect_equal(data$get_num_bins(), 3)
  expect_equal(data$get_num_sequences(TRUE), 5)
  expect_equal(data$get_num_sequences(), 111)
  expect_equal(data$get_num_samples(), 0)
  expect_equal(data$get_num_treatments(), 0)

  table <- export_dataset(data)
  dataset2 <- import_dataset(table)

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

  table <- export_dataset(just_bins)

  dataset <- import_dataset(table)

  expect_equal(dataset$get_num_bins("otu"), 531)
  expect_equal(dataset$get_num_sequences(), 113963)

  expect_equal(dataset$get_dataset_name(), just_bins$get_dataset_name())

  expect_equal(dataset$get_num_sequences(), just_bins$get_num_sequences())
  expect_equal(dataset$get_num_treatments(), just_bins$get_num_treatments())
  expect_equal(dataset$get_num_samples(), just_bins$get_num_samples())
  expect_equal(length(get_sequence_names(dataset)), 0)

  expect_error(import_dataset(table, c("sequence_data")))
  expect_error(import_dataset(table, c("metadata")))
  expect_error(import_dataset(table, c("references")))
  expect_error(import_dataset(table, c("sequence_tree")))

  attr(table, "rdataset_version") <- "0.2.1"
  expect_error(import_dataset(table))
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

  table <- export_dataset(just_seqs)

  dataset <- import_dataset(table)

  expect_equal(dataset$get_num_bins("otu"), 0)
  expect_equal(dataset$get_num_sequences(), 113963)
  expect_equal(just_seqs$get_num_sequences(TRUE), 2425)

  expect_equal(dataset$get_dataset_name(), just_seqs$get_dataset_name())

  expect_equal(dataset$get_num_sequences(), just_seqs$get_num_sequences())
  expect_equal(dataset$get_num_treatments(), just_seqs$get_num_treatments())
  expect_equal(dataset$get_num_samples(), just_seqs$get_num_samples())
  expect_equal(length(get_sequence_names(dataset)), 2425)

  expect_error(import_dataset(table, c("bin_data")))
  expect_error(import_dataset(table, c("sample_tree")))
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

  table <- export_dataset(just_bins)

  data <- import_dataset(table)

  expect_equal(data$get_num_bins("otu"), 531)
  expect_equal(data$get_num_sequences(), 113963)

  expect_equal(data$get_dataset_name(), just_bins$get_dataset_name())

  expect_equal(data$get_num_sequences(), just_bins$get_num_sequences())
  expect_equal(data$get_num_treatments(), just_bins$get_num_treatments())
  expect_equal(data$get_num_samples(), just_bins$get_num_samples())
  expect_equal(length(get_sequence_names(data)), 0)

  data <- dataset$new()

  table <- export_dataset(data, c("sequence_data", "bin_data"))
  expect_equal(length(table), 0)
})

test_that("import - with tags", {
  miseq <- miseq_sop_example()

  expect_equal(miseq$get_num_bins("otu"), 531)
  expect_equal(miseq$get_num_sequences(), 113963)
  expect_equal(miseq$get_num_sequences(TRUE), 2425)
  expect_equal(nrow(get_bin_representative_sequences(miseq)), 531)

  # just export bin data, no sequence data
  table <- export_dataset(miseq)

  just_bins <- import_dataset(table, c("bin_data"))

  expect_equal(just_bins$get_num_bins("otu"), 531)
  expect_equal(just_bins$get_num_sequences(), 113963)
  expect_equal(just_bins$get_num_sequences(TRUE), 0)

  expect_equal(get_dataset_name(just_bins), get_dataset_name(miseq))
  expect_equal(just_bins$get_num_treatments(), 0)
  expect_equal(just_bins$get_num_samples(), 0)
  expect_equal(length(get_sequence_names(just_bins)), 0)
  expect_equal(nrow(get_bin_representative_sequences(just_bins)), 0)
})
