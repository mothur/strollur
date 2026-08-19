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

test_that("test report - type = sequence_tree and sample_tree", {
  tree <- ape::read.tree(strollur_example("final.phylip.tre.gz"))

  data <- new_dataset()
  add(data, table = tree, type = "sequence_tree")

  tree <- report(data, type = "sequence_tree")

  expect_equal(sort(names(data, "sequence")), sort(tree$tip.label))
  expect_equal(tree$edge[1:5, 1], c(2426, 2427, 2427, 2426, 2428))
  expect_equal(tree$edge[1:5, 2], c(2427, 1, 2, 2428, 2429))
  expect_equal(
    round(tree$edge.length[1:5], digits = 3),
    c(NaN, 0.004, 0.004, 0.000, 0.002)
  )

  dataset_t <- read_mothur(
    fasta = strollur_example("final.fasta.gz"),
    count = strollur_example("final.count_table.gz"),
    otu_list = strollur_example("final.opti_mcc.list.gz"),
    dataset_name = "miseq_sop"
  )

  sample_tree <- ape::read.tree(
    strollur_example("final.opti_mcc.jclass.ave.tre")
  )

  add(dataset_t, table = sample_tree, type = "sample_tree")

  tree <- report(dataset_t, type = "sample_tree")

  expect_equal(sort(names(dataset_t, "sample")), sort(tree$tip.label))
  expect_equal(tree$edge[1:5, 1], c(20, 21, 22, 23, 24))

  xdev_remove_samples(dataset_t, c("F3D1", "F3D141"))

  tree <- report(dataset_t, type = "sample_tree")

  expect_equal(sort(names(dataset_t, "sample")), sort(tree$tip.label))
})
