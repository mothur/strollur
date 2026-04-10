# tests add function

test_that("test add - errors", {
  data <- new_dataset()

  # not a strollur object
  x <- 10
  expect_error(add(x))

  fasta_data <- read_fasta(strollur_example("final.fasta.gz"))

  # test bad type
  expect_error(add(data, table = fasta_data, type = bad_type))

  # test bad report_type
  expect_error(
    add(data, table = fasta_data, type = reports),
    "'report_type' is required when adding a report."
  )

  add(data, table = fasta_data)
  expect_equal(count(data, type = "sequences"), 2425)
})
