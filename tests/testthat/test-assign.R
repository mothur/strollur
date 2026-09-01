# tests assign function

test_that("test assign - piping", {
  data <- new_dataset()

  table <- data.frame(
    sequence_name = c("seq1", "seq2", "seq3"),
    abundance = c(40, 50, 10)
  )

  abunds <- assign(data, table, type = sequence_abundance) |> abundance()

  expect_equal(abunds$sequence_name, c("seq1", "seq2", "seq3"))
  expect_equal(abunds$abundance, c(40, 50, 10))
})

test_that("test assign - errors", {
  data <- new_dataset()

  # not a strollur object
  x <- 10
  expect_error(assign(x))

  table <- readRDS(strollur_example("miseq_abundance_by_sample.rds"))

  # test no quotes on type
  assign(data, table = table, type = sequence_abundance)

  expect_equal(count(data, type = "sequence"), 113963)

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

  expect_equal(count(data, type = "bin"), 531)

  assign(data,
    table = bin_reps,
    type = "bin_representative", bin_type = "otu"
  )

  expect_equal(nrow(report(data, type = bin_representative)), 531)
})

test_that("test assign - sample_distance", {
  # Create a new empty strollur object named 'example_dataset'
  data <- strollur::new_dataset(dataset_name = "example_dataset")

  shared_file <- strollur_example("final.opti_mcc.shared")
  dist_file <- strollur_example("final.opti_mcc.jclass.0.03.column.dist")
  reference <- strollur::new_reference(name = "jclass estimator distances",
                                       vendor = "mothur_v1.48.6")

  df <- read_mothur_shared(shared_file)
  xdev_assign_bins(data, table = df, bin_type = "otu")

  sample_dists <- readr::read_table(dist_file,
                                    col_names = FALSE,
                                    show_col_types = FALSE)
  strollur::assign(data, table = sample_dists, type = "sample_distance")

  expect_equal(count(data, type = "sample"), 19)

  dists <- xdev_get_sample_distances(data)

  expect_equal(nrow(dists), 171)
  expect_equal(ncol(dists), 3)
  expect_equal(round(dists[1:3, 3], 2), c(0.47, 0.53, 0.55))

})
