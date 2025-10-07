test_that("write_mothur_taxonomy - errors", {
  expect_error(write_mothur_taxonomy("Bad_type"))

  # no file name with nameless dataset
  data <- dataset$new()
  expect_error(write_mothur_taxonomy(data))
})

test_that("write_mothur_taxonomy", {
  miseq <- miseq_sop_example()

  output <- write_mothur_taxonomy(miseq, get_full_name("miseq.taxonomy"))

  df <- read_mothur_taxonomy(output)

  remove_file(output)

  expect_equal(df[[1]], miseq$get_sequence_names())
  tax1 <- paste0(
    "Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);",
    "Lachnospiraceae(85);Lachnospiraceae_unclassified(85);"
  )

  tax2 <- paste0(
    "Bacteria(100);\"Bacteroidetes\"(100);\"Bacteroidia\"(99);",
    "\"Bacteroidales\"(99);\"Porphyromonadaceae\"(90);",
    "\"Porphyromonadaceae\"_unclassified(90);"
  )

  tax3 <- paste0(
    "Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);",
    "Lachnospiraceae(99);Lachnospiraceae_unclassified(99);"
  )

  expect_equal(df[[2]][1], tax1)
  expect_equal(df[[2]][100], tax2)
  expect_equal(df[[2]][50], tax3)

  data <- read_mothur(
    otu_list = rdataset_example("final.opti_mcc.list"),
    dataset_name = "data"
  )

  expect_equal(write_mothur_taxonomy(data), "no_sequence_taxonomy")
})
