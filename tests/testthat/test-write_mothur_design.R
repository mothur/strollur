test_that("write_mothur_design - errors", {
  expect_error(write_mothur_design("Bad_type"))

  # no file name with nameless dataset
  data <- dataset$new()
  expect_error(write_mothur_design(data))
})

test_that("write_mothur_design", {
  miseq <- miseq_sop_example()

  output <- write_mothur_design(miseq, get_full_name("miseq.design"))

  data <- read_mothur(
    otu_shared = rdataset_example("final.opti_mcc.shared"),
    design = output
  )

  remove_file(output)

  expect_equal(count(data, "treatments"), 2)
  expect_equal(count(data, "samples"), 19)

  data <- read_mothur(
    otu_list = rdataset_example("final.opti_mcc.list"),
    dataset_name = "data"
  )

  expect_equal(write_mothur_design(data), "no_design_data")
})
