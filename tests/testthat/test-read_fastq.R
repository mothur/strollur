# test read_fastq

test_that("test read_fastq - errors", {
  expect_error(read_fastq("non_existant_filename"))
})

test_that("test read_fastq", {
  fastq_data <- strollur::read_fastq("tiny.fastq.gz", strollur_example())

  names <- c(
    "M00967:43:000000000-A3JHG:1:1101:18327:1699",
    "M00967:43:000000000-A3JHG:1:1101:14069:1827",
    "M00967:43:000000000-A3JHG:1:1101:18044:1900"
  )

  expect_equal(fastq_data$sequence_name, names)

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

  expect_equal(fastq_data$sequence, sequences)

  score_1_15 <- c(2, 29, 29, 32, 32, 33, 32, 33, 33, 37, 37, 37, 38, 38, 38)
  score_2_15 <- c(18, 32, 32, 30, 32, 33, 33, 35, 33, 37, 37, 33, 36, 38, 38)
  score_3_15 <- c(33, 32, 31, 33, 33, 33, 32, 33, 33, 37, 37, 37, 38, 38, 38)

  expect_equal(fastq_data$quality_score[[1]][1:15], score_1_15)
  expect_equal(fastq_data$quality_score[[2]][1:15], score_2_15)
  expect_equal(fastq_data$quality_score[[3]][1:15], score_3_15)

  fastq_data <- strollur::read_fastq("tiny.fastq.gz",
                                     strollur_example(),
                                     format = "solexa")

  score_1_15 <- c(0, 2, 2, 3, 3, 4, 3, 4, 4, 6, 6, 6, 7, 7, 7)
  score_2_15 <- c(0, 3, 3, 2, 3, 4, 4, 5, 4, 6, 6, 4, 5, 7, 7)
  score_3_15 <- c(4, 3, 3, 4, 4, 4, 3, 4, 4, 6, 6, 6, 7, 7, 7)

  expect_equal(fastq_data$quality_score[[1]][1:15], score_1_15)
  expect_equal(fastq_data$quality_score[[2]][1:15], score_2_15)
  expect_equal(fastq_data$quality_score[[3]][1:15], score_3_15)

  # negative scores
  expect_error(strollur::read_fastq("tiny.fastq.gz",
                                     strollur_example(),
                                     format = "illumina"))

  fastq_data <- strollur::read_fastq("tiny.fastq.gz",
                                    strollur_example(),
                                    format = "sanger")

  score_1_15 <- c(2, 29, 29, 32, 32, 33, 32, 33, 33, 37, 37, 37, 38, 38, 38)
  score_2_15 <- c(18, 32, 32, 30, 32, 33, 33, 35, 33, 37, 37, 33, 36, 38, 38)
  score_3_15 <- c(33, 32, 31, 33, 33, 33, 32, 33, 33, 37, 37, 37, 38, 38, 38)

  expect_equal(fastq_data$quality_score[[1]][1:15], score_1_15)
  expect_equal(fastq_data$quality_score[[2]][1:15], score_2_15)
  expect_equal(fastq_data$quality_score[[3]][1:15], score_3_15)
})

test_that("test convert_qual", {
  quality_string <- "3AA?ABBDBF"
  scores <- convert_qual(quality_string)

  expect_equal(scores, c(18, 32, 32, 30, 32, 33, 33, 35, 33, 37))

  quality_string <- "FBEGGEGGGG"
  scores <- convert_qual(quality_string)
  expect_equal(scores, c(37, 33, 36, 38, 38, 36, 38, 38, 38, 38))

  quality_string <- "HGGCFDEFGGGH"
  scores <- convert_qual(quality_string, format = "solexa")
  expect_equal(scores, c(8, 7, 7, 4, 6, 5, 5, 6, 7, 7, 7, 8))
})
