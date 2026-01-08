test_that("write_mothur_count - errors", {
  expect_error(write_mothur_count("Bad_type"))

  # no file name with nameless dataset
  data <- dataset$new()
  expect_error(write_mothur_count(data))
})

test_that("write_mothur_count - with sample data", {
  miseq <- miseq_sop_example()

  output <- write_mothur_count(miseq, get_full_name("miseq.taxonomy"))

  data <- read_mothur(count = output)

  remove_file(output)

  expected <- abundance(miseq, type = "sequences", by_sample = TRUE)

  # remove treatment column since data does not include treatments
  expected <- expected[, -c(4)]

  actual <- abundance(data, type = "sequences", by_sample = TRUE)

  # sort to same order
  sorted_indices <- order(expected[[1]])
  expected <- expected[sorted_indices, ]
  sorted_indices <- order(actual[[1]])
  actual <- actual[sorted_indices, ]

  # remove row names since sort made then not match
  rownames(actual) <- NULL
  rownames(expected) <- NULL

  expect_equal(expected, actual)

  data <- read_mothur(
    otu_shared = rdataset_example("final.opti_mcc.shared"),
    dataset_name = "data"
  )

  expect_equal(write_mothur_count(data), "no_sequence_data")
})

test_that("write_mothur_count - without sample data", {
  data2 <- read_mothur(count = rdataset_example("test_nogroups.count_table"))

  output <- write_mothur_count(data2, get_full_name("data2.taxonomy"))

  data <- read_mothur(count = output)

  remove_file(output)

  expected <- abundance(data2, , type = "sequences", by_sample = TRUE)
  actual <- abundance(data, type = "sequences", by_sample = TRUE)

  # sort to same order
  sorted_indices <- order(expected[[1]])
  expected <- expected[sorted_indices, ]
  sorted_indices <- order(actual[[1]])
  actual <- actual[sorted_indices, ]

  # remove row names since sort made then not match
  rownames(actual) <- NULL
  rownames(expected) <- NULL

  expect_equal(expected, actual)
})
