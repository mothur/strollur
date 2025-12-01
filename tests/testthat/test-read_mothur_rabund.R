# test read_mothur_rabund

test_that("test read_mothur_rabund - errors", {
  expect_error(read_mothur_rabund("non_existant_filename"))
})

test_that("test read_mothur_rabund", {
  table <- read_mothur_rabund(rdataset_example("final.opti_mcc.rabund"))

  data <- dataset$new()
  assign_bins(data, table)

  expect_equal(num(data), 2425)
  expect_equal(num(data, "bins", "otu"), 531)
})
