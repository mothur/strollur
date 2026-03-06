# tests count function

test_that("count - non_dataset", {
  # test non dataset class with names
  df <- data.frame(var = c(1, 3), tmp = c(2, 4))
  expect_equal(count(df)$n, 2)
})

test_that("count - samples & treatments", {
  miseq <- miseq_sop_example()

  expect_equal(count(miseq, type = "samples"), 19)
  expect_equal(count(miseq, type = "treatments"), 2)

  expect_error(count(miseq, type = "not_valid"))
})

test_that("count - bins", {
  miseq <- miseq_sop_example()

  expect_equal(count(miseq, type = "bins", bin_type = "otu"), 531)
  expect_equal(count(miseq, type = "bins", bin_type = "asv"), 2425)
  expect_equal(count(miseq, type = "bins", bin_type = "phylotype"), 63)

  # with samples parameter - numbins with seqs from both "F3D0" && "F3D1"
  expect_equal(count(miseq, "bins", samples = c("F3D0", "F3D1")), 125)

  # with samples parameter - distinct = true
  # numbins exclusive to both "F3D0" && "F3D1"
  expect_equal(count(miseq, "bins", "otu", c("F3D0", "F3D1"), TRUE), 1)

  # bad sample
  message <- capture_output(count(
    miseq,
    type = "bins",
    samples = c("bad_sample", "F3D1")
  ))
  expect_true(grepl("samples requested, ignoring.", message))
})

test_that("count - sequences", {
  miseq <- miseq_sop_example()

  # all seqs
  expect_equal(count(miseq, type = "sequences"), 113963)
  # unique seqs
  expect_equal(count(miseq, type = "sequences", distinct = TRUE), 2425)

  # with samples parameter - distinct = false
  expect_equal(count(miseq, samples = c("F3D0", "F3D1")), 9385)

  # with samples parameter - distinct = true
  expect_equal(
    count(miseq, "sequences", "w", c("F3D0", "F3D1"), TRUE),
    2
  )

  data <- read_mothur(otu_shared = strollur_example("final.opti_mcc.shared"))
  expect_equal(count(data, samples = c("F3D0", "F3D1")), 10351)
})
