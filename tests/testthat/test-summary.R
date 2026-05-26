# tests summary function

test_that("summary tests", {
  # test non dataset type
  # Should we just error for non dataset types?
  # contigs_report <- readRDS(strollur_example("miseq_contigs_report.rds"))
  # expect_equal(class(summary(contigs_report)), "table")

  miseq <- miseq_sop_example()

  # To get the summary of your FASTA data
  df <- summary(data = miseq, type = "sequence", verbose = FALSE)
  headers <- c(
    "starts", "ends", "nbases", "ambigs",
    "polymers", "numns", "numseqs"
  )
  expect_true(all(headers %in% names(df)))

  starts <- rep(1, 8)
  expect_equal(starts, df$starts)

  ends <- rep(375, 8)
  expect_equal(ends, df$ends)
  expect_equal(df$nbases[[7]], 256)

  ambigs <- rep(0, 8)
  expect_equal(ambigs, df$ambigs)
  # numNs
  expect_equal(ambigs, df$numns)

  expect_equal(df$polymers[[4]], 4)
  expect_equal(df$numseqs[[4]], 56982.5)

  # summarize contigs_report
  df <- summary(
    data = miseq, type = "report",
    report_type = "contigs_report", verbose = FALSE
  )

  headers <- c(
    "Expected_Errors", "Length", "MisMatches", "Num_Ns",
    "Overlap_End", "Overlap_Length", "Overlap_Start"
  )
  expect_true(all(headers %in% names(df)))
  expect_equal(df$MisMatches[[7]], 120)
  expect_equal(df$Overlap_Length[[1]], 232)
  expect_equal(df$MisMatches[[4]], 1)
  expect_equal(df$Overlap_End[[6]], 253)
  expect_equal(df$Num_Ns[[3]], 0)
  expect_equal(df$Length[[2]], 252)

  # remove sample 'F3D0' to produce a scrap report
  xdev_remove_samples(data = miseq, samples = c("F3D0"))

  # summarize FASTA data after removal of sample F3D0
  df <- summary(data = miseq, type = "sequence")

  expect_equal(starts, df[[1]])
  expect_equal(ends, df[[2]])
  expect_equal(df[[7]][4], 53887)

  # summarize scrapped data -
  # sequences and bins scrapped by removing the sample "F3D0"
  df <- summary(data = miseq, type = "scrap", verbose = FALSE)

  expect_equal(names(df), c("type", "trash_code", "unique", "total"))
  expect_equal(df[[1]], c("sequence", "otu", "asv", "phylotype"))
  expect_equal(df[[2]], rep("remove_samples", 4))
  expect_equal(df[[3]], c(101, 14, 101, 2))
  expect_equal(df[[4]], c(109, 14, 109, 2))
})

test_that("add, assign, abundance error tests", {
  contigs_report <- readRDS(strollur_example("miseq_contigs_report.rds"))

  # no report type provided
  expect_error(add(
    data = data, table = contigs_report, type = "report",
    table_names = list(sequence_name = "Name")
  ))

  # invalid type
  expect_error(add(
    data = data, table = contigs_report, type = "bin",
    table_names = list(sequence_name = "Name")
  ))

  # invalid type
  expect_error(assign(
    data = data, table = contigs_report, type = "report",
    table_names = list(sequence_name = "Name")
  ))

  # invalid type
  expect_error(abundance(
    data = data, table = contigs_report,
    type = "report"
  ))
})
