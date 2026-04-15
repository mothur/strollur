# test read_mothur_rabund

test_that("test read_mothur_rabund - errors", {
  expect_error(read_mothur_rabund("non_existant_filename"))
})

test_that("test read_mothur_rabund", {
  table <- read_mothur_rabund(strollur_example("final.opti_mcc.rabund"))

  data <- new_dataset()
  assign(data = data, table = table, type = "bins")

  expect_equal(count(data), 2425)
  expect_equal(count(data, "bins", "otu"), 531)
})
