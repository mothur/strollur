# adds coverage for internal (xint) and developer (xdev) Rcpp functions

test_that("xdev_abundance", {
  data <- new_dataset()

  # bad sample
  message <- capture_output(xdev_abundance(data, type = "bad_type"))
  expect_true(grepl(
    "bad_type is not a valid type for the abundance function",
    message
  ))
})

test_that("xdev_assign_bin_taxonomy", {
  data <- new_dataset()

  # can't assign bin taxonomy when you have not assigned bins
  expect_error(xdev_assign_bin_taxonomy(data, data.frame()))
})

test_that("xdev_assign_sequence_abundance", {
  data <- new_dataset()

  add(data, data.frame(sequence_names = c("seq1")))

  table <- data.frame(
    sequence_names = c("seq2", "seq3"),
    abundances = c(100, 200)
  )

  # must assign abundances for all sequences
  expect_error(xdev_assign_sequence_abundance(data, table))
})

test_that("xdev_get_by_sample", {
  data <- new_dataset()

  # can't assign bin taxonomy when you have not assigned bins
  expect_error(xdev_get_by_sample(data, "badType"))
})

xdev_get_by_sample
test_that("xdev_assign_bins, assign with reference", {
  data <- new_dataset()

  # bad sample
  reference <- new_reference("myReference")

  table <- data.frame(
    bin_names = c("bin1", "bin2", "bin2"),
    sequence_names = c("seq1", "seq2", "seq3")
  )

  xdev_assign_bins(data, table, reference = reference)

  ref_report <- report(data, type = "references")

  expect_equal(nrow(ref_report), 1)

  xdev_assign_bin_representative_sequences(
    data,
    data.frame(
      bin_names = c("bin1", "bin2"),
      sequence_names = c("seq1", "seq2")
    ),
    reference = reference
  )

  ref_report <- report(data, type = "references")

  expect_equal(nrow(ref_report), 2)
})
