# tests import of dataset object

test_that("import - miseq_sop_example", {
  # full dataset
  miseq <- miseq_sop_example()

  expect_equal(count(miseq, "bins", "phylotype"), 63)
  expect_equal(count(miseq, distinct = TRUE), 2425)

  phylo_bins_to_remove <- c("Phylo01", "Phylo02")
  reasons_to_remove <- c("testing", "testing")

  # remove some bins to allow for filtering
  xdev_remove_bins(
    miseq, phylo_bins_to_remove,
    reasons_to_remove, "phylotype"
  )

  expect_equal(count(miseq, "bins", bin_type = "phylotype"), 61)
  expect_equal(count(miseq, distinct = TRUE), 825)
  expect_equal(count(miseq), 39177)

  exported_miseq <- export_dataset(miseq)

  expect_equal(sum(exported_miseq$sequence_data$include_sequence), 825)

  dataset_t <- import_dataset(exported_miseq)

  expect_equal(names(dataset_t, "dataset"), names(miseq, "dataset"))
  expect_equal(count(dataset_t, distinct = TRUE), count(miseq, distinct = TRUE))
  expect_equal(count(dataset_t), count(miseq))
  expect_equal(count(dataset_t, "treatments"), count(miseq, "treatments"))
  expect_equal(count(dataset_t, "samples"), count(miseq, "samples"))
  expect_equal(count(dataset_t, "bins", "otu"), count(miseq, "bins", "otu"))
  expect_equal(
    count(dataset_t, "bins", "phylotype"), count(miseq, "bins", "phylotype")
  )
  expect_equal(count(dataset_t, "bins", "asv"), count(miseq, "bins", "asv"))
  expect_equal(
    totals(dataset_t),
    totals(miseq)
  )
  expect_equal(
    totals(dataset_t, "treatments"),
    totals(miseq, "treatments")
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

  expect_equal(count(data, "bins"), 3)
  expect_equal(count(data, distinct = TRUE), 5)
  expect_equal(count(data), 111)
  expect_equal(count(data, "samples"), 0)
  expect_equal(count(data, "treatments"), 0)

  table <- export_dataset(data)
  dataset2 <- import_dataset(table)

  expect_equal(count(dataset2, "bins"), 3)
  expect_equal(count(dataset2, distinct = TRUE), 5)
  expect_equal(count(dataset2), 111)
  expect_equal(count(dataset2, "samples"), 0)
  expect_equal(count(dataset2, "treatments"), 0)
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

  expect_equal(count(just_bins, "bins"), 531)
  expect_equal(count(just_bins), 113963)

  table <- export_dataset(just_bins)

  dataset <- import_dataset(table)

  expect_equal(count(dataset, "bins"), 531)
  expect_equal(count(dataset), 113963)

  expect_equal(
    names(dataset, "dataset"),
    names(just_bins, "dataset")
  )

  expect_equal(count(dataset), count(just_bins))
  expect_equal(count(dataset, "treatments"), count(just_bins, "treatments"))
  expect_equal(count(dataset, "samples"), count(just_bins, "samples"))
  expect_equal(length(names(dataset, "sequences")), 0)

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

  expect_equal(count(just_seqs, "bins"), 0)
  expect_equal(count(just_seqs), 113963)
  expect_equal(count(just_seqs, distinct = TRUE), 2425)

  table <- export_dataset(just_seqs)

  dataset <- import_dataset(table)

  expect_equal(count(dataset, "bins"), 0)
  expect_equal(count(dataset), 113963)
  expect_equal(count(dataset, distinct = TRUE), 2425)

  expect_equal(names(dataset, "dataset"), names(just_seqs, "dataset"))

  expect_equal(count(dataset), count(just_seqs))
  expect_equal(count(dataset, "treatments"), count(just_seqs, "treatments"))
  expect_equal(count(dataset, "samples"), count(just_seqs, "samples"))
  expect_equal(length(names(dataset)), 2425)

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

  expect_equal(count(just_bins, "bins", "otu"), 531)
  expect_equal(count(just_bins), 113963)

  table <- export_dataset(just_bins)

  data <- import_dataset(table)

  expect_equal(count(data, "bins", "otu"), 531)
  expect_equal(count(data), 113963)

  expect_equal(names(data, "dataset"), names(just_bins, "dataset"))

  expect_equal(count(data), count(just_bins))
  expect_equal(count(data, "treatments"), count(just_bins, "treatments"))
  expect_equal(count(data, "samples"), count(just_bins, "samples"))
  expect_equal(length(names(data)), 0)

  data <- dataset$new()

  table <- export_dataset(data, c("sequence_data", "bin_data"))
  expect_equal(length(table), 0)
})

test_that("import - with tags", {
  miseq <- miseq_sop_example()

  expect_equal(count(miseq, "bins", "otu"), 531)
  expect_equal(count(miseq), 113963)
  expect_equal(count(miseq, distinct = TRUE), 2425)
  expect_equal(nrow(get_bin_representative_sequences(miseq)), 531)

  # just export bin data, no sequence data
  table <- export_dataset(miseq)

  just_bins <- import_dataset(table, c("bin_data"))

  expect_equal(count(just_bins, "bins", "otu"), 531)
  expect_equal(count(just_bins), 113963)
  expect_equal(count(just_bins, distinct = TRUE), 0)

  expect_equal(names(miseq, "dataset"), names(just_bins, "dataset"))
  expect_equal(count(just_bins, "treatments"), 0)
  expect_equal(count(just_bins, "samples"), 0)
  expect_equal(length(names(just_bins)), 0)
  expect_equal(nrow(get_bin_representative_sequences(just_bins)), 0)
})
