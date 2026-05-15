# test read_mothur_shared

test_that("test read_mothur_shared - errors", {
  expect_error(read_mothur_shared("non_existant_filename"))
})

test_that("test read_mothur_shared", {
  shared_data <- read_mothur_shared(strollur_example("final.opti_mcc.shared"))

  data <- new_dataset()
  assign(data = data, table = shared_data, type = "bin")

  expect_equal(count(data), 113963)
  expect_equal(count(data, "sample"), 19)
  expect_equal(count(data, "bin", "otu"), 531)
})
