# test read_fastq

test_that("test read_quality - errors", {
  expect_error(read_quality("non_existant_filename"))
})

test_that("test read_quality", {
  qual_data <- strollur::read_quality("tiny.qual", path = strollur_example())

  names <- c(
    "M00967:43:000000000-A3JHG:1:1101:18327:1699",
    "M00967:43:000000000-A3JHG:1:1101:14069:1827",
    "M00967:43:000000000-A3JHG:1:1101:18044:1900"
  )

  expect_equal(qual_data$sequence_name, names)

  score_1_15 <- c(2, 29, 29, 32, 32, 33, 32, 33, 33, 37, 37, 37, 38, 38, 38)
  score_2_15 <- c(18, 32, 32, 30, 32, 33, 33, 35, 33, 37, 37, 33, 36, 38, 38)
  score_3_15 <- c(33, 32, 31, 33, 33, 33, 32, 33, 33, 37, 37, 37, 38, 38, 38)

  expect_equal(qual_data$quality_score[[1]][1:15], score_1_15)
  expect_equal(qual_data$quality_score[[2]][1:15], score_2_15)
  expect_equal(qual_data$quality_score[[3]][1:15], score_3_15)
})
