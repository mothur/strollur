# test "write_mothur_rabund"

test_that("write_mothur_rabund - errors", {
  expect_error(write_mothur_rabund("Bad_type"))

  # no file name with nameless dataset
  data <- new_dataset()
  expect_error(write_mothur_rabund(data))
})

test_that("write_mothur_rabund", {
  miseq <- miseq_sop_example()
  bin_types <- miseq$get_bin_types()

  file_root <- get_full_name("test-miseq")

  outputs <- c(
    get_full_name("test-miseq.otu.rabund"),
    get_full_name("test-miseq.asv.rabund"),
    get_full_name("test-miseq.phylotype.rabund")
  )

  expect_equal(write_mothur_rabund(miseq, file_root), outputs)

  df <- read_mothur_rabund(outputs[1])

  expected <- xdev_abundance(
    data = miseq, type = "bin",
    bin_type = bin_types[1]
  )
  names(expected) <- c("bin_name", "abundance")

  expect_equal(ncol(df), ncol(expected))
  expect_equal(df$abundance, expected$abundance)

  # cleanup
  for (output in outputs) {
    remove_file(output)
  }
})
