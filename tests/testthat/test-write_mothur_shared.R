# test "write_mothur_shared"

test_that("write_mothur_shared - errors", {
  expect_error(write_mothur_shared("Bad_type"))

  # no file name with nameless dataset
  data <- new_dataset()
  expect_error(write_mothur_shared(data))
})

test_that("write_mothur_shared", {
  miseq <- miseq_sop_example()
  bin_types <- miseq$get_bin_types()

  file_root <- get_full_name("test-miseq")

  outputs <- c(
    "test-miseq.otu.shared",
    "test-miseq.asv.shared",
    "test-miseq.phylotype.shared"
  )
  outputs <- paste0(normalizePath(test_path()), .Platform$file.sep, outputs)

  expect_equal(write_mothur_shared(miseq, file_root), outputs)

  df <- read_mothur_shared(outputs[1])
  names(df) <- c("bin_names", "abundances", "samples")

  expected <- xdev_abundance(
    data = miseq, type = "bins",
    bin_type = bin_types[1], by_sample = TRUE
  )

  # remove treatment column
  expected <- expected[, -c(4)]
  names(expected) <- c("bin_names", "abundances", "samples")

  expect_equal(df, expected)

  # cleanup
  for (output in outputs) {
    remove_file(output)
  }
})
