# test read_mothur_shared

test_that("test read_mothur_shared - errors", {
  expect_error(read_mothur_shared("non_existant_filename"))
})

test_that("test read_mothur_shared", {
  shared_data <- read_mothur_shared(rdataset_example("final.opti_mcc.shared"))

  data <- dataset$new()
  assign_bins(data, shared_data)

  expect_equal(num(data), 113963)
  expect_equal(num(data, "samples"), 19)
  expect_equal(num(data, "bins", "otu"), 531)
})
