# tests add function

test_that("test add - piping", {
  data <- new_dataset()

  table <- data.frame(sequence_name = c("seq1", "seq2", "seq3"))

  abunds <- add(data, table) |> abundance()

  expect_equal(abunds$sequence_name, c("seq1", "seq2", "seq3"))
  expect_equal(abunds$abundance, rep(1, 3))
})

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
    add(data, table = fasta_data, type = report),
    "'report_type' is required when adding a report."
  )

  add(data, table = fasta_data)
  expect_equal(count(data, type = "sequence"), 2425)

  reference <- readr::read_csv(strollur_example("references.csv"),
    col_names = TRUE, show_col_types = FALSE
  )

  add(data, reference, "resource_reference")

  references <- report(data, "resource_reference")

  # random spot checks
  expect_equal(nrow(references), 2)
  expect_equal(references[[1, "name"]], "R phylotypr package")
  expect_equal(references[[2, "name"]], "silva.bacteria.fasta")
  expect_equal(references[[1, "note"]], "classification using Bayesian method")
  expect_equal(
    references[[2, "note"]],
    "alignment reference trimmed to V4 region"
  )
})
