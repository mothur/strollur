# tests abundance function

test_that("abundance", {
  # not a strollur object
  x <- 10
  expect_error(abundance(x))

  # without ""

  strollur <- new_dataset()
  expect_equal(nrow(abundance(strollur, type = bins, bin_type = asv)), 0)
})
