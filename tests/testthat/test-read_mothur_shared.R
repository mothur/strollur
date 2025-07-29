# test read_mothur_shared

test_that("test read_mothur_shared - errors", {
  expect_error(read_mothur_shared("non_existant_filename"))
})

test_that("test read_mothur_shared", {
  shared_data <- read_mothur_shared(rdataset_example("final.opti_mcc.shared"))

  dataset <- sequence_data$new()
  dataset$assign_bins(
    shared_data$bin_id,
    shared_data$abundance,
    shared_data$sample
  )

  expect_equal(dataset$get_num_sequences(), 113963)
  expect_equal(dataset$get_num_samples(), 19)
  expect_equal(dataset$get_num_bins("otu"), 531)
})
