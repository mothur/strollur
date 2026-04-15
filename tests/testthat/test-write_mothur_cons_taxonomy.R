# test "write_mothur_cons_taxonomy"

test_that("write_mothur_cons_taxonomy - errors", {
  expect_error(write_mothur_cons_taxonomy("Bad_type"))

  # no file name with nameless dataset
  data <- new_dataset()
  expect_error(write_mothur_cons_taxonomy(data))
})

test_that("write_mothur_cons_taxonomy", {
  miseq <- miseq_sop_example()
  bin_types <- miseq$get_bin_types()

  file_root <- get_full_name("test-miseq")

  outputs <- c(
    "test-miseq.otu.cons.taxonomy",
    "test-miseq.asv.cons.taxonomy",
    "test-miseq.phylotype.cons.taxonomy"
  )
  outputs <- paste0(normalizePath(test_path()), .Platform$file.sep, outputs)

  expect_equal(outputs, write_mothur_cons_taxonomy(miseq, file_root))

  df <- read_mothur_cons_taxonomy(outputs[1])

  expect_equal(df[[1]], xdev_names(miseq, "bins", bin_types[1]))
  expect_equal(df[[2]], xdev_abundance(
    data = miseq, type = "bins",
    bin_type = bin_types[1]
  )[[2]])
  tax1 <- paste0(
    "Bacteria(100);\"Bacteroidetes\"(100);\"Bacteroidia\"(100);",
    "\"Bacteroidales\"(100);\"Porphyromonadaceae\"(100);",
    "\"Porphyromonadaceae\"_unclassified(100);"
  )

  expect_equal(df[[3]][1], tax1)

  df <- read_mothur_cons_taxonomy(outputs[2])

  expect_equal(df[[1]], xdev_names(miseq, "bins", bin_types[2]))
  expect_equal(df[[2]], xdev_abundance(
    data = miseq, type = "bins",
    bin_type = bin_types[2]
  )[[2]])
  tax2 <- paste0(
    "Bacteria(100);Firmicutes(100);Clostridia(100);",
    "Clostridiales(100);Lachnospiraceae(100);",
    "Lachnospiraceae_unclassified(100);"
  )

  expect_equal(df[[3]][100], tax2)

  df <- read_mothur_cons_taxonomy(outputs[3])

  expect_equal(df[[1]], xdev_names(miseq, "bins", bin_types[3]))
  expect_equal(df[[2]], xdev_abundance(
    data = miseq, type = "bins",
    bin_type = bin_types[3]
  )[[2]])
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
