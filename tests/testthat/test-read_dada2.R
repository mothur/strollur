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

  expect_equal(names(data, "dataset"), "data")
  expect_equal(names(data, "samples"), c("sample1", "sample2", "sample3"))
  expect_equal(count(data, "sequences", distinct = TRUE), 100)
  expect_equal(count(data, "sequences"), 45150)
  expect_equal(get_sequences(data), rep("ATGCTTT", 100))
  expect_equal(length(names(data, "sequences")), 100)
})
