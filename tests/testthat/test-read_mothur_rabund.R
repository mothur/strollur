# test read_mothur_rabund

test_that("test read_mothur_rabund - errors", {
  expect_error(read_mothur_rabund("non_existant_filename"))
})

test_that("test read_mothur_rabund", {
  data <- read_mothur_rabund(rdataset_example("final.opti_mcc.rabund"))

  dataset <- sequence_data$new()
  dataset$assign_bins(data)

  expect_equal(dataset$get_num_sequences(), 2425)
  expect_equal(dataset$get_num_bins("otu"), 531)
})
