test_that("write_fastq - errors", {
  expect_error(write_fastq("Bad_type"))

  # no file name with nameless dataset
  data <- new_dataset()
  expect_error(write_fastq(data))
})

test_that("write_fastq", {
  data <- strollur::new_dataset(dataset_name = "example_dataset")

  fastq_data <-
    strollur::read_fastq(fastq = strollur_example("tiny.fastq.gz"))

  strollur::add(data = data, table = fastq_data, type = "fastq")

  output <- write_fastq(data, get_full_name("example_dataset"))

  df <- read_fastq(output)

  remove_file(output)

  expect_equal(df[[1]], names(data, "sequence"))

  seq_1 <- paste0(
    "NACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGC",
    "CTGCCAAGTCAGCGGTAAAATTGCGGGGCTCAACCCCGTACAGCCGTTGAAACTGC",
    "CGGGCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATA",
    "GATATCACGCAGAACCCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTG",
    "AGGCACGAAAGTGCGGGGATCAAACAG"
  )

  seq_2 <- paste0(
    "TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGC",
    "CTGCCAAGTCAGCGGTAAAATTGCGGGGCTCAACCCCGTACAGCCGTTGAAACTGC",
    "CGGGCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATA",
    "GATATCACGCAGAACCCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTG",
    "AGGCACGAAAGTGCGGGGATCAAACAG"
  )

  seq_3 <- paste0(
    "TACGGAGGATGCGAGCGTTGTCCGGAATCACTGGGCGTAAAGGGCGCGTAGGCGGTT",
    "TAATAAGTCAGTGGTGAAAACTGAGGGCTCAACCCTCAGCCTGCCACTGATACTGTT",
    "AGACTTGAGTATGGAAGAGGAGAATGGAATTCCTAGTGTAGCGGTGAAATGCGTAGA",
    "TATTAGGAGGAACACCAGTGGCGAAGGCGATTCTCTGGGCCAAGACTGACGCTGAGG",
    "CGCGAAAGCGTGGGGAGCAAACA"
  )
  sequences <- c(seq_1, seq_2, seq_3)

  expect_equal(df$sequence, sequences)

  score_1_15 <- c(2, 29, 29, 32, 32, 33, 32, 33, 33, 37, 37, 37, 38, 38, 38)
  score_2_15 <- c(18, 32, 32, 30, 32, 33, 33, 35, 33, 37, 37, 33, 36, 38, 38)
  score_3_15 <- c(33, 32, 31, 33, 33, 33, 32, 33, 33, 37, 37, 37, 38, 38, 38)

  expect_equal(df$quality_score[[1]][1:15], score_1_15)
  expect_equal(df$quality_score[[2]][1:15], score_2_15)
  expect_equal(df$quality_score[[3]][1:15], score_3_15)

  output <- write_fastq(data, get_full_name("example_dataset"))

  df <- read_fastq(output, format = "solexa")

  remove_file(output)

  score_1_15 <- c(0, 2, 2, 3, 3, 4, 3, 4, 4, 6, 6, 6, 7, 7, 7)
  score_2_15 <- c(0, 3, 3, 2, 3, 4, 4, 5, 4, 6, 6, 4, 5, 7, 7)
  score_3_15 <- c(4, 3, 3, 4, 4, 4, 3, 4, 4, 6, 6, 6, 7, 7, 7)

  expect_equal(df$quality_score[[1]][1:15], score_1_15)
  expect_equal(df$quality_score[[2]][1:15], score_2_15)
  expect_equal(df$quality_score[[3]][1:15], score_3_15)

  data <- read_mothur(
      otu_list = strollur_example("final.opti_mcc.list.gz"),
      dataset_name = "data"
  )

  expect_equal(write_fastq(data), "no_fastq_data")
})
