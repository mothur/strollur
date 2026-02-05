# test read_mothur_count

test_that("test read_mothur_count compressed with samples", {
  # read compressed count file with samples
  results <- read_mothur_count(strollur_example("test.count_table"))

  expect_equal(c("sequence_names", "samples", "abundances"), colnames(results))

  abunds <- rep(1, 1000)
  expect_equal(abunds, results$abundance)

  first_ten <- c(
    "M00967_43_000000000-A3JHG_1_1101_10386_25574",
    "M00967_43_000000000-A3JHG_1_1101_10579_16062",
    "M00967_43_000000000-A3JHG_1_1101_12161_19944",
    "M00967_43_000000000-A3JHG_1_1101_12704_3254",
    "M00967_43_000000000-A3JHG_1_1101_12778_23425",
    "M00967_43_000000000-A3JHG_1_1101_12954_6506",
    "M00967_43_000000000-A3JHG_1_1101_13090_11008",
    "M00967_43_000000000-A3JHG_1_1101_13440_22261",
    "M00967_43_000000000-A3JHG_1_1101_13874_18116",
    "M00967_43_000000000-A3JHG_1_1101_15924_15103"
  )

  expect_equal(first_ten, results$sequence_names[1:10])

  first_ten_samples <- c(
    "F3D2", "F3D146", "F3D150", "F3D145",
    "F3D147", "F3D8", "F3D147", "F3D150",
    "F3D150", "F3D148"
  )

  expect_equal(first_ten_samples, results$samples[1:10])
})

test_that("test read_mothur_count uncompressed with samples", {
  # read uncompressed count file with samples
  results <- read_mothur_count(strollur_example("test.full.count_table"))

  expect_equal(c("sequence_names", "samples", "abundances"), colnames(results))

  abunds <- rep(1, 1000)
  expect_equal(abunds, results$abundances)

  first_ten <- c(
    "M00967_43_000000000-A3JHG_1_1101_10386_25574",
    "M00967_43_000000000-A3JHG_1_1101_10579_16062",
    "M00967_43_000000000-A3JHG_1_1101_12161_19944",
    "M00967_43_000000000-A3JHG_1_1101_12704_3254",
    "M00967_43_000000000-A3JHG_1_1101_12778_23425",
    "M00967_43_000000000-A3JHG_1_1101_12954_6506",
    "M00967_43_000000000-A3JHG_1_1101_13090_11008",
    "M00967_43_000000000-A3JHG_1_1101_13440_22261",
    "M00967_43_000000000-A3JHG_1_1101_13874_18116",
    "M00967_43_000000000-A3JHG_1_1101_15924_15103"
  )

  expect_equal(first_ten, results$sequence_names[1:10])

  first_ten_samples <- c(
    "F3D2", "F3D146", "F3D150", "F3D145",
    "F3D147", "F3D8", "F3D147", "F3D150",
    "F3D150", "F3D148"
  )

  expect_equal(first_ten_samples, results$samples[1:10])
})

test_that("test read_mothur_count uncompressed NO samples", {
  # read uncompressed count file without samples
  results <- read_mothur_count(
    strollur_example("test_nogroups.count_table")
  )

  expect_equal(c("sequence_names", "abundances"), colnames(results))

  abunds <- c(4773, 14378, 108939, 1734, 150122, 15750)
  expect_equal(abunds, results$abundance)

  names <- c(
    "chimera1", "non_chimera2ParentB", "chimera1ParentB",
    "non_chimera2", "chimera1ParentA", "non_chimera2ParentA"
  )

  expect_equal(names, results$sequence_names)
})
