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
    abundance(dataset_t, type = "samples"),
    abundance(miseq, type = "samples")
  )
  expect_equal(
    abundance(dataset_t, "treatments"),
    abundance(miseq, "treatments")
  )

  dfd <- report(dataset_t, "bin_representatives")
  dfm <- report(miseq, "bin_representatives")

  for (otu in dfm[[1]]) {
    expect_equal(
      dfd |> filter(otu_names == otu),
      dfm |> filter(otu_names == otu)
    )
  }

  data <- strollur$new()
  abunds <- c(1, 10, 100)
  bins <- c("otu1", "otu2", "otu3")

  xdev_assign_bins(
    data = data,
    table = data.frame(bin_names = bins, abundances = abunds)
  )

  expect_equal(count(data, "bins"), 3)
  expect_equal(count(data, distinct = TRUE), 3)
  expect_equal(count(data), 111)
  expect_equal(count(data, "samples"), 0)
  expect_equal(count(data, "treatments"), 0)

  table <- export_dataset(data)
  dataset2 <- import_dataset(table)

  expect_equal(count(dataset2, "bins"), 3)
  expect_equal(count(dataset2, distinct = TRUE), 3)
  expect_equal(count(dataset2), 111)
  expect_equal(count(dataset2, "samples"), 0)
  expect_equal(count(dataset2, "treatments"), 0)
})

test_that("import - no sequence data", {
  # just shared and constax dataset
  just_bins <- read_mothur(
    otu_shared = strollur_example("final.opti_mcc.shared"),
    cons_taxonomy = strollur_example("final.cons.taxonomy"),
    design = strollur_example("mouse.time.design"),
    sample_tree = strollur_example("final.opti_mcc.jclass.ave.tre"),
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
  expect_equal(length(names(dataset, "sequences")), 531)

  expect_error(import_dataset(table, c("sequence_data")))
  expect_error(import_dataset(table, c("metadata")))
  expect_error(import_dataset(table, c("references")))
  expect_error(import_dataset(table, c("sequence_tree")))

  attr(table, "strollur_version") <- "0.2.1"
  expect_error(import_dataset(table))
})

test_that("import - no bin data", {
  # just shared and constax dataset
  just_seqs <- read_mothur(
    fasta = strollur_example("final.fasta.gz"),
    count = strollur_example("final.count_table.gz"),
    taxonomy = strollur_example("final.taxonomy.gz"),
    design = strollur_example("mouse.time.design"),
    sequence_tree = strollur_example("final.phylip.tre.gz"),
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
    otu_shared = strollur_example("final.opti_mcc.shared"),
    cons_taxonomy = strollur_example("final.cons.taxonomy"),
    design = strollur_example("mouse.time.design"),
    sample_tree = strollur_example("final.opti_mcc.jclass.ave.tre"),
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
  expect_equal(length(names(data)), 531)

  data <- strollur$new()

  table <- export_dataset(data)
  expect_equal(length(table), 0)
})
