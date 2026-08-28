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

test_that("test assign - quality", {
  # Create a new empty strollur object named 'example_dataset'
  data <- strollur::new_dataset(dataset_name = "example_dataset")

  # Read quality data into data.frame
  quality_data <-
    strollur::read_quality(quality = strollur_example("tiny.qual"))
  fasta_data <-
    strollur::read_fasta(fasta = strollur_example("tiny.fasta"))

  # Add FASTA sequence data
  strollur::add(data = data, table = fasta_data, type = "sequence")
  strollur::assign(data, table = quality_data, type = "quality")

  quality_report <- strollur::report(data, type = "quality")

  names <- c(
    "M00967:43:000000000-A3JHG:1:1101:18327:1699",
    "M00967:43:000000000-A3JHG:1:1101:14069:1827",
    "M00967:43:000000000-A3JHG:1:1101:18044:1900"
  )

  expect_equal(quality_report$sequence_name, names)

  score_1_15 <- c(2, 29, 29, 32, 32, 33, 32, 33, 33, 37, 37, 37, 38, 38, 38)
  score_2_15 <- c(18, 32, 32, 30, 32, 33, 33, 35, 33, 37, 37, 33, 36, 38, 38)
  score_3_15 <- c(33, 32, 31, 33, 33, 33, 32, 33, 33, 37, 37, 37, 38, 38, 38)

  expect_equal(quality_report$quality_score[[1]][1:15], score_1_15)
  expect_equal(quality_report$quality_score[[2]][1:15], score_2_15)
  expect_equal(quality_report$quality_score[[3]][1:15], score_3_15)
})
