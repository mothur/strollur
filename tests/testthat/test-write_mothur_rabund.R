# test "write_mothur_rabund"

test_that("write_mothur_rabund - errors", {
  expect_error(write_mothur_rabund("Bad_type"))

  # no file name with nameless dataset
  data <- dataset$new()
  expect_error(write_mothur_rabund(data))
})

test_that("write_mothur_rabund", {
  miseq <- miseq_sop_example()
  bin_types <- miseq$get_bin_types()

  file_root <- get_full_name("test-miseq")

  outputs <- c(
    "test-miseq.otu.rabund",
    "test-miseq.asv.rabund",
    "test-miseq.phylotype.rabund"
  )
  outputs <- paste0(normalizePath(test_path()), .Platform$file.sep, outputs)

  expect_equal(write_mothur_rabund(miseq, file_root), outputs)

  df <- read_mothur_rabund(outputs[1])

  expected <- abundance(data = miseq, type = "bins", bin_type = bin_types[1])
  names(expected) <- c("bin_names", "abundances")

  expect_equal(ncol(df), ncol(expected))
  expect_equal(df$abundances, expected$abundances)

  # cleanup
  for (output in outputs) {
    remove_file(output)
  }
})
