# test "write_mothur_cons_taxonomy"

test_that("write_mothur_cons_taxonomy - errors", {
  expect_error(write_mothur_cons_taxonomy("Bad_type"))

  # no file name with nameless dataset
  data <- dataset$new()
  expect_error(write_mothur_cons_taxonomy(data))

  # no bin classifications
  data <- read_mothur(
    otu_list = rdataset_example("final.opti_mcc.list"),
    dataset_name = "temp"
  )
})

test_that("write_mothur_cons_taxonomy", {
  miseq <- miseq_sop_example()
  bin_types <- miseq$get_bin_types()

  outputs <- write_mothur_cons_taxonomy(miseq)

  df <- read_mothur_cons_taxonomy(outputs[1])

  expect_equal(df[[1]], miseq$get_bin_names(bin_types[1]))
  expect_equal(df[[2]], get_rabund_vector(miseq$data, bin_types[1]))
  tax1 <- paste0(
    "Bacteria(100);\"Bacteroidetes\"(100);\"Bacteroidia\"(100);",
    "\"Bacteroidales\"(100);\"Porphyromonadaceae\"(100);",
    "\"Porphyromonadaceae\"_unclassified(100);"
  )

  expect_equal(df[[3]][1], tax1)

  df <- read_mothur_cons_taxonomy(outputs[2])

  expect_equal(df[[1]], miseq$get_bin_names(bin_types[2]))
  expect_equal(df[[2]], get_rabund_vector(miseq$data, bin_types[2]))
  tax2 <- paste0(
    "Bacteria(100);Firmicutes(100);Clostridia(100);",
    "Clostridiales(100);Lachnospiraceae(100);",
    "Lachnospiraceae_unclassified(100);"
  )

  expect_equal(df[[3]][100], tax2)

  df <- read_mothur_cons_taxonomy(outputs[3])

  expect_equal(df[[1]], miseq$get_bin_names(bin_types[3]))
  expect_equal(df[[2]], get_rabund_vector(miseq$data, bin_types[3]))
  tax3 <- paste0(
    "Bacteria(100);Firmicutes(100);Clostridia(100);",
    "Clostridiales(100);Lachnospiraceae(100);",
    "Marvinbryantia(100);"
  )

  expect_equal(df[[3]][50], tax3)

  # cleanup
  for (output in outputs) {
    remove_file(output)
  }
})
