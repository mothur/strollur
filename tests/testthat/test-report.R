# tests report function

test_that("test report - errors", {
  data <- new_dataset()

  # add references, custom report and metadata
  contigs_report <- readr::read_tsv(
    strollur_example(
      "final.contigs_report.gz"
    ),
    col_names = TRUE, show_col_types = FALSE
  )
  add(data,
    table = contigs_report, type = "reports",
    report_type = "contigs_report", list(sequence_name = "Name")
  )

  # bad sample
  message <- capture_output(report(
    data,
    type = "bad_type"
  ))
  expect_true(grepl("does not include a report named", message))
})

test_that("test report - type = fasta", {
  data <- new_dataset()

  names <- c("seq1", "seq2", "seq3")
  seqs <- c("..ATGCTTGC..", ".AT-GC-TT-GC.", "ATGCTTGC")
  comments <- c("sdf", "wer", "tyu")

  add(data, table = data.frame(
    sequence_names = names,
    sequences = seqs,
    comments = comments
  ))

  fasta <- report(data, type = "fasta")

  expect_equal(ncol(fasta), 3)
  expect_equal(fasta$sequence_names, names)
  expect_equal(fasta$sequences, seqs)
  expect_equal(fasta$comments, comments)

  data <- new_dataset()

  add(data, table = data.frame(
    sequence_names = names,
    sequences = seqs
  ))

  fasta <- report(data, type = "fasta")

  expect_equal(ncol(fasta), 2)
  expect_equal(fasta$sequence_names, names)
  expect_equal(fasta$sequences, seqs)
})
