# test read_mothur_shared

test_that("test read_mothur_shared - errors", {
  expect_error(read_mothur_shared("non_existant_filename"))
})

test_that("test read_mothur_shared", {
  shared_data <- read_mothur_shared(rdataset_example("final.opti_mcc.shared"))

  data <- dataset$new()
  assign_bins(data, shared_data)

  expect_equal(data$get_num_sequences(), 113963)
  expect_equal(data$get_num_samples(), 19)
  expect_equal(data$get_num_bins("otu"), 531)
})
