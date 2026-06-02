# tests report function

test_that("test report - errors", {
  data <- new_dataset()

  # add references, custom report and metadata
  contigs_report <- readRDS(strollur_example("miseq_contigs_report.rds"))
  add(data,
    table = contigs_report, type = "report",
    report_type = "contigs_report", list(sequence_name = "Name")
  )

  # bad sample
  message <- capture_output(report(
    data,
    type = "bad_type"
  ))
  expect_true(grepl("does not include a report named", message))

  # not a strollur object
  x <- 10
  expect_error(report(x))
})

test_that("test report - type = fasta", {
  data <- new_dataset()

  names <- c("seq1", "seq2", "seq3")
  seqs <- c("..ATGCTTGC..", ".AT-GC-TT-GC.", "ATGCTTGC")
  comments <- c("sdf", "wer", "tyu")

  add(data, table = data.frame(
    sequence_name = names,
    sequence = seqs,
    comment = comments
  ))

  fasta <- report(data, type = "fasta")

  expect_equal(ncol(fasta), 3)
  expect_equal(fasta$sequence_name, names)
  expect_equal(fasta$sequence, seqs)
  expect_equal(fasta$comment, comments)

  xdev_remove_sequences(data,
    sequence_names = c("seq1"),
    trash_tags = c("scrap_report_test")
  )

  sequence_scrap_report <- report(data, "sequence_scrap")

  expect_equal(ncol(sequence_scrap_report), 2)
  expect_equal(names(sequence_scrap_report), c("sequence_name", "trash_code"))
  expect_equal(sequence_scrap_report$sequence_name, c("seq1"))
  expect_equal(sequence_scrap_report$trash_code, c("scrap_report_test"))

  data <- new_dataset()

  sequence_scrap_report <- report(data, "sequence_scrap")
  expect_equal(sequence_scrap_report, data.frame())

  add(data, table = data.frame(
    sequence_name = names,
    sequence = seqs
  ))

  fasta <- report(data, type = "fasta")

  expect_equal(ncol(fasta), 2)
  expect_equal(fasta$sequence_name, names)
  expect_equal(fasta$sequence, seqs)
})
