# tests names function

test_that("names - dataset", {
  data <- new_dataset("soil")

  expect_equal(names(data, type = "dataset"), "soil")

  # test non dataset class with names
  df <- data.frame(var = c(1, 3), tmp = c(2, 4))

  expect_equal(names(df), c("var", "tmp"))

  expect_error(names(data, type = "not_valid_type"))
})

test_that("names - samples & treatments & reports", {
  miseq <- miseq_sop_example()

  samples <- c(
    "F3D0", "F3D1", "F3D141", "F3D142", "F3D143",
    "F3D144", "F3D145", "F3D146", "F3D147", "F3D148", "F3D149",
    "F3D150", "F3D2", "F3D3", "F3D5",
    "F3D6", "F3D7", "F3D8", "F3D9"
  )

  expect_equal(names(miseq, type = "sample"), samples)

  treatments <- c("Early", "Late")

  expect_equal(names(miseq, type = "treatment"), treatments)

  reports <- c("contigs_report")

  expect_equal(names(miseq, type = "report"), reports)
})

test_that("names - bins", {
  miseq <- miseq_sop_example()

  expect_equal(length(names(miseq, type = "bin", bin_type = "otu")), 531)
  expect_equal(length(names(miseq, type = "bin", bin_type = "asv")), 2425)
  expect_equal(
    length(names(miseq, type = "bin", bin_type = "phylotype")),
    63
  )

  # with samples parameter
  expect_equal(length(names(miseq, "bin", samples = c("F3D0", "F3D1"))), 125)

  # with samples parameter - distinct = true
  expect_equal(length(names(miseq, "bin", "otu", c("F3D0", "F3D1"), TRUE)), 1)

  # bad sample
  message <- capture_output(names(
    miseq,
    type = "bin",
    samples = c("bad_sample", "F3D1")
  ))
  expect_true(grepl("samples requested, ignoring.", message))
})

test_that("names - sequences", {
  miseq <- miseq_sop_example()

  expect_equal(length(names(miseq, type = "sequence")), 2425)

  # with samples parameter - distinct = false
  expect_equal(length(names(miseq, samples = c("F3D0", "F3D1"))), 133)

  # with samples parameter - distinct = true
  expect_equal(
    length(names(miseq, "sequence", "w", c("F3D0", "F3D1"), TRUE)),
    2
  )
})
