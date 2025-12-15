# test read_mothur_shared

test_that("test read_mothur_shared - errors", {
  expect_error(read_mothur_shared("non_existant_filename"))
})

test_that("test read_mothur_shared", {
  shared_data <- read_mothur_shared(rdataset_example("final.opti_mcc.shared"))

  data <- dataset$new()
  assign(data = data, table = shared_data, type = "bins")

  expect_equal(count(data), 113963)
  expect_equal(count(data, "samples"), 19)
  expect_equal(count(data, "bins", "otu"), 531)
})
