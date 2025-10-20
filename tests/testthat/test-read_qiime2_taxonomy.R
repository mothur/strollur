# test read_qiime2_taxonomy

test_that("test read_qiime2_taxonomy - errors", {
  expect_error(read_qiime2_taxonomy("non_existant_filename"))
  expect_error(read_qiime2_taxonomy(
    rdataset_example("table.qza"),
    file_root <- get_full_name("test-qiime2")
  ))
  unlink(get_full_name("test-qiime2"), recursive = TRUE)
})

test_that("test read_qiime2_taxonomy", {
  tax_data <- read_qiime2_taxonomy(rdataset_example("taxonomy.qza"))

  expect_equal(nrow(tax_data$data), 759)
  expect_equal(ncol(tax_data$data), 3)

  tax <- paste0(
    "k__Bacteria; p__Proteobacteria; c__Betaproteobacteria;",
    " o__Burkholderiales"
  )
  expect_equal(tax_data$data[264, 2], tax)
  expect_equal(tax_data$data[264, 1], "3ef461f213bfc675ddf0b1bfeef2dc52")
  expect_equal(round(tax_data$data[264, 3] * 100), 97)
})
