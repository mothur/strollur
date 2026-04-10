# test "write_mothur_list"

test_that("write_mothur_list - errors", {
  expect_error(write_mothur_list("Bad_type"))

  # no file name with nameless dataset
  data <- strollur$new()
  expect_error(write_mothur_list(data))
})

test_that("write_mothur_list", {
  miseq <- miseq_sop_example()
  bin_types <- miseq$get_bin_types()

  file_root <- get_full_name("test-miseq")

  outputs <- c(
    "test-miseq.otu.list",
    "test-miseq.asv.list",
    "test-miseq.phylotype.list"
  )
  outputs <- paste0(normalizePath(test_path()), .Platform$file.sep, outputs)

  expect_equal(write_mothur_list(miseq, file_root), outputs)

  df <- read_mothur_list(outputs[1])
  names(df) <- c("otu_id", "seq_id")

  expect_equal(df, xdev_report(miseq, "sequence_bin_assignments", bin_types[1]))

  # cleanup
  for (output in outputs) {
    remove_file(output)
  }
})
