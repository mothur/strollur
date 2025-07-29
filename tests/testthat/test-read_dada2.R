# test read_dada2

test_that("test read_dada2 - errors", {
  expect_error(read_dada2())
})

test_that("test read_dada2", {
  # 100 seqs, 3 samples
  dada2_seqtab <- matrix(
    c(1:300),
    nrow = 3,
    ncol = 100,
    byrow = TRUE
  )

  rownames(dada2_seqtab) <- c("sample1", "sample2", "sample3")
  colnames(dada2_seqtab) <- rep("ATGCTTT", 100)

  data <- read_dada2(dada2_seqtab, "data")

  expect_equal(data$data$dataset_name, "data")
  expect_equal(data$get_samples(), c("sample1", "sample2", "sample3"))
  expect_equal(data$get_num_sequences(TRUE), 100)
  expect_equal(data$get_num_sequences(), 45150)
  expect_equal(data$get_sequences(), rep("ATGCTTT", 100))
  expect_equal(length(data$get_ids()), 100)
})
