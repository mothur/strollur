# tests assign function

test_that("test assign - errors", {
  data <- new_dataset()

  # not a strollur object
  x <- 10
  expect_error(assign(x))

  table <- readRDS(strollur_example("miseq_abundance_by_sample.rds"))

  # test no quotes on type
  assign(data, table = table, type = sequence_abundance)

  expect_equal(count(data, type = "sequences"), 113963)

  expect_error(assign(data, table = table, type = "bad_type"))
})

test_that("test assign - bin reps", {
  bin_reps <- readRDS(strollur_example("miseq_representative_sequences.rds"))
  otu_data <- readRDS(strollur_example("miseq_shared_otu.rds"))

  # create dataset
  data <- read_mothur(
    fasta = strollur_example("final.fasta.gz"),
    count = strollur_example("final.count_table.gz"),
    otu_list = strollur_example("final.opti_mcc.list.gz"),
    dataset_name = "test"
  )

  expect_equal(count(data, type = "bins"), 531)

  assign(data,
    table = bin_reps,
    type = "bin_representatives", bin_type = "otu"
  )

  expect_equal(nrow(report(data, type = bin_representatives)), 531)
})
